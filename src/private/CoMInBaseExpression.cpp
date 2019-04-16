/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/CoMInBaseExpression.h>
#include <unordered_set>

namespace DynamicalPlanner {
    namespace Private {
        class ComInBaseHessianEvaluable;
        class CoMInBaseJacobianEvaluable;
        class CoMInBasePositionEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::ComInBaseHessianEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    std::unordered_set<size_t> m_nonZeros;

public:

    ComInBaseHessianEvaluable(ExpressionsServer* es, long column)
        : levi::DefaultEvaluable(3, es->jointsPosition().rows(), "d^2(CoM_B)/d(q^2)[:," + std::to_string(column) + "]")
    {

        assert(es);
        std::string baseFrame = es->getFloatingBase();
        const iDynTree::Model& model = es->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(es->getFloatingBase());
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex baseLink = model.getFrameLink(baseFrameIndex);
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
        assert(ok);

        double totalMass = 0.0;

        for(size_t l=0; l < model.getNrOfLinks(); l++)
        {
            totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
        }

        m_cols.resize(model.getNrOfDOFs(), levi::Null(3,1));

        iDynTree::LinkConstPtr linkPtr;
        iDynTree::IJointConstPtr jointPtr;
        iDynTree::LinkIndex visitedLink, childLink, parentLink;
        size_t jointIndex;

        iDynTree::LinkIndex columnChild = traversal.getChildLinkIndexFromJointIndex(model, static_cast<iDynTree::JointIndex>(column));
        std::string columnChildName = model.getLinkName(columnChild);
        iDynTree::LinkIndex columnParent = traversal.getParentLinkIndexFromJointIndex(model, static_cast<iDynTree::JointIndex>(column));

        levi::Expression thisS = es->motionSubSpaceAsCrossProduct(static_cast<iDynTree::JointIndex>(column), columnParent, columnChild);

        std::unordered_set<size_t> otherJointsInSameBranch;
        std::string linkName;

        for (iDynTree::LinkIndex l = 0; l < static_cast<int>(model.getNrOfLinks()); ++l) {

            assert(model.isValidLinkIndex(l));
            linkPtr = model.getLink(l);
            linkName = model.getLinkName(l);
            visitedLink = l;
            levi::Scalar relativeMass(linkPtr->getInertia().getMass() / totalMass, " m_" + linkName + "/M");
            levi::Constant linkComPosition(iDynTree::toEigen(linkPtr->getInertia().getCenterOfMass()), "p_CoM_" + linkName);
            levi::Expression linkComPosition6D = levi::Expression::Vertcat(linkComPosition, levi::Null(3,1), "[p_CoM_" + linkName + ";0]");

            std::vector<size_t> visitedJoints;

            while(visitedLink != baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                if (jointIndex == static_cast<size_t>(column)) {

                    levi::Expression jacobianColumnTopRows = es->adjointTransform(linkName, columnChildName).block(0,0,3,6) * es->motionSubSpaceVector(jointPtr->getIndex(), columnParent, columnChild);
                    levi::Expression jacobianColumn = levi::Expression::Vertcat(jacobianColumnTopRows,
                                                                                levi::Null(3,1), "[" + jacobianColumnTopRows.name() + ";0]");

                    levi::Expression constantPart = es->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), columnParent, columnChild) *
                            es->adjointTransform(model.getLinkName(columnChild), linkName) * linkComPosition6D +
                        es->adjointTransform(model.getLinkName(columnChild), linkName) * jacobianColumn;

                    while (visitedLink != baseLink) {
                        jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                        jointIndex = static_cast<size_t>(jointPtr->getIndex());

                        childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
                        parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                        m_cols[jointIndex] = m_cols[jointIndex] +
                            relativeMass * es->adjointTransform(baseFrame, model.getLinkName(childLink)).block(0,0,3,6) *
                                es->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), parentLink, childLink) *
                                es->adjointTransform(model.getLinkName(childLink), columnChildName)* constantPart;

                        m_nonZeros.insert(jointIndex);

                        visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
                    }

                    otherJointsInSameBranch.insert(visitedJoints.begin(), visitedJoints.end());

                } else {
                    visitedJoints.push_back(jointIndex);
                    visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
                }
            }
        }

        m_nonZeros.insert(otherJointsInSameBranch.begin(), otherJointsInSameBranch.end());

        for (size_t otherCols : otherJointsInSameBranch) {
            m_cols[otherCols] = es->comInBaseHessian(static_cast<long>(otherCols)).col(column); //here we are exploiting the fact that d_sj*d_sk is equal to
                                                                                                //d_sk*d_sj. Thus we compute only the derivatives for which
                                                                                                //s_k is closer to the base than s_j
        }

        addDependencies(es->jointsPosition());

    }

//    virtual levi::ColumnExpression col(Eigen::Index col) final {
//        return m_cols[static_cast<size_t>(col)];
//    }

//    virtual levi::ScalarExpression element(Eigen::Index row, Eigen::Index col) final {
//        return m_cols[static_cast<size_t>(col)](row, 0);
//    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        for (size_t j : m_nonZeros) {
            m_evaluationBuffer.col(static_cast<Eigen::Index>(j)) = m_cols[j].evaluate();
        }

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        for (auto& expression : m_cols) {
            expression.clearDerivativesCache();
        }
    }

};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
ComInBaseHessianEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    return m_cols[static_cast<size_t>(column)].getColumnDerivative(0, variable);
}

class DynamicalPlanner::Private::CoMInBaseJacobianEvaluable : public levi::DefaultEvaluable {
    ExpressionsServer* m_expressionsServer;
    iDynTree::MatrixDynSize m_jacobian;
    levi::Variable m_jointsVariable;


//    std::vector<levi::Expression> m_columns;

public:

    CoMInBaseJacobianEvaluable(ExpressionsServer* expressionsServer)
        : levi::DefaultEvaluable(3, expressionsServer->jointsPosition().rows(), "B_J_[B, CoM]")
          , m_expressionsServer(expressionsServer)
          , m_jointsVariable(expressionsServer->jointsPosition())
    {
        m_jacobian.resize(3, static_cast<unsigned int>(expressionsServer->jointsPosition().rows()));
        m_jacobian.zero();
//        const iDynTree::Model& model = expressionsServer->model();

//        iDynTree::Traversal traversal;
//        std::string baseFrame = expressionsServer->getFloatingBase();
//        iDynTree::LinkIndex baseLink = model.getLinkIndex(baseFrame);
//        assert(baseLink != iDynTree::LINK_INVALID_INDEX);
//        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
//        assert(ok);

//        double totalMass = 0.0;

//        for(size_t l=0; l < model.getNrOfLinks(); l++)
//        {
//            totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
//        }

//        m_columns.resize(static_cast<size_t>(expressionsServer->jointsPosition().rows()), levi::Null(3,1));

//        iDynTree::IJointConstPtr jointPtr;
//        iDynTree::LinkIndex linkIndex;
//        iDynTree::LinkConstPtr linkPtr;
//        size_t jointIndex;
//        std::string linkName;

//        iDynTree::LinkIndex visitedLink, childLink, parentLink;

//        Eigen::Matrix<double, 6,1> motionSubSpaceVector;

//        for (size_t l = 0; l < model.getNrOfLinks(); ++l) {

//            linkIndex = static_cast<iDynTree::LinkIndex>(l);
//            assert(model.isValidLinkIndex(linkIndex));
//            linkName = model.getLinkName(linkIndex);
//            linkPtr = model.getLink(linkIndex);
//            levi::Scalar relativeMass(linkPtr->getInertia().getMass() / totalMass, " m_" + linkName + "/M");
//            levi::Constant linkComPosition(iDynTree::toEigen(linkPtr->getInertia().getCenterOfMass()), "p_CoM_" + linkName);
//            levi::Expression linkComPosition6D = levi::Expression::Vertcat(linkComPosition, levi::Null(3,1), "[p_CoM_" + linkName + ";0]");

//            visitedLink = linkIndex;

//            while (visitedLink != baseLink) {
//                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
//                jointIndex = static_cast<size_t>(jointPtr->getIndex());

//                childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
//                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

//                levi::Expression jacobianColumnTopRows = expressionsServer->adjointTransform(linkName, model.getLinkName(childLink)).block(0,0,3,6) * expressionsServer->motionSubSpaceVector(jointPtr->getIndex(), parentLink, childLink);
//                levi::Expression jacobianColumn = levi::Expression::Vertcat(jacobianColumnTopRows,
//                                                                            levi::Null(3,1), "[" + jacobianColumnTopRows.name() + ";0]");

//                m_columns[jointIndex] = m_columns[jointIndex] +
//                    relativeMass * expressionsServer->adjointTransform(baseFrame, model.getLinkName(childLink)).block(0,0,3,6) * (expressionsServer->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), parentLink, childLink) *
//                                                                                                                                     expressionsServer->adjointTransform(model.getLinkName(childLink), linkName) * linkComPosition6D +
//                                                                                                                                     expressionsServer->adjointTransform(model.getLinkName(childLink), linkName) * jacobianColumn);
//                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
//            }
//        }

        addDependencies(m_jointsVariable);
    }

//    virtual levi::ColumnExpression column(Eigen::Index col) final {
//        return m_columns[static_cast<size_t>(col)];
//    }

//    virtual levi::ScalarExpression element(Eigen::Index row, Eigen::Index col) final {
//        return m_columns[static_cast<size_t>(col)](row, 0);
//    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        SharedKinDynComputationsPointer kinDyn = m_expressionsServer->currentKinDyn();

        iDynTree::Vector4 normalizedQuaternion;
        iDynTree::toEigen(normalizedQuaternion) = iDynTree::toEigen(m_expressionsServer->currentState().base_quaternion).normalized();
        iDynTree::Rotation baseRotation;
        baseRotation.fromQuaternion(normalizedQuaternion);

        bool ok = kinDyn->getCenterOfMassJacobian(m_expressionsServer->currentState(), m_jacobian);
        assert(ok);

        m_evaluationBuffer = iDynTree::toEigen(baseRotation.inverse()) * iDynTree::toEigen(m_jacobian).rightCols(m_expressionsServer->currentState().s.size());

        return m_evaluationBuffer;
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
//        for (auto& expression : m_columns) {
//            expression.clearDerivativesCache();
//        }
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
CoMInBaseJacobianEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows()) { //m_expression contains the jointsVariable specified in the constructor

        return m_expressionsServer->comInBaseHessian(column);

    } else {
        return levi::Null(6, variable->dimension());
    }

}

class DynamicalPlanner::Private::CoMInBasePositionEvaluable : public levi::DefaultEvaluable {

    ExpressionsServer* m_expressionsServer;
    levi::Expression m_derivative;
    levi::Variable m_jointsVariable;

public:

    CoMInBasePositionEvaluable(ExpressionsServer* expressionsServer)
        : levi::DefaultEvaluable( 3, 1, "b_p_CoM")
          , m_expressionsServer(expressionsServer)
          , m_jointsVariable(expressionsServer->jointsPosition())
    {

        m_derivative = levi::ExpressionComponent<CoMInBaseJacobianEvaluable>(expressionsServer);
        addDependencies(m_jointsVariable);

    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        SharedKinDynComputationsPointer kinDyn = m_expressionsServer->currentKinDyn();

        iDynTree::Vector4 normalizedQuaternion;
        iDynTree::toEigen(normalizedQuaternion) = iDynTree::toEigen(m_expressionsServer->currentState().base_quaternion).normalized();
        iDynTree::Rotation baseRotation;
        baseRotation.fromQuaternion(normalizedQuaternion);

        iDynTree::Transform baseTransform;
        baseTransform.setRotation(baseRotation);
        iDynTree::Position basePosition;
        iDynTree::toEigen(basePosition) = iDynTree::toEigen(m_expressionsServer->currentState().base_position);
        baseTransform.setPosition(basePosition);

        iDynTree::Position comPosition = kinDyn->getCenterOfMassPosition(m_expressionsServer->currentState());

        m_evaluationBuffer = iDynTree::toEigen(baseTransform.inverse() * comPosition);

        return m_evaluationBuffer;
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_derivative.clearDerivativesCache();
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
CoMInBasePositionEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows()) { //m_expression contains the jointsVariable specified in the constructor

        assert(column == 0);

        levi::unused(column);

        return m_derivative;
    } else {
        return levi::Null(3, variable->dimension());
    }

}

levi::Expression DynamicalPlanner::Private::CoMInBaseExpression(ExpressionsServer *expressionsServer)
{
    return levi::ExpressionComponent<CoMInBasePositionEvaluable>(expressionsServer);
}

levi::Expression DynamicalPlanner::Private::CoMInBaseJointsDoubleDerivative(ExpressionsServer *expressionsServer, long column)
{
    return levi::ExpressionComponent<ComInBaseHessianEvaluable>(expressionsServer, column);
}

levi::Expression DynamicalPlanner::Private::CoMAdjointTransformWrench(const levi::Expression &worldToBaseRotation, const levi::Expression &comInBasePosition)
{
    levi::Null topRightBlock(3,3);
    levi::Expression leftCols = levi::Expression::Vertcat(worldToBaseRotation, worldToBaseRotation * (-1.0 * comInBasePosition).skew(), "leftCols");
    levi::Expression rightCols = levi::Expression::Vertcat(topRightBlock, worldToBaseRotation, "rightCols");
    return levi::Expression::Horzcat(leftCols, rightCols, "Gbar_X*_B");
}
