/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/MomentumInBaseExpression.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <unordered_set>
#include <cassert>


namespace DynamicalPlanner {
    namespace Private {
        class MomentumInBaseEvaluable;
        class MomentumInBaseJointsDerivativeEvaluable;
        class MomentumInBaseBaseTwistDerivativeEvaluable;
        class MomentumInBaseJointsDoubleDerivativeEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

/*-------------------------- Momentum In Base -----------------------------------------*/

class DynamicalPlanner::Private::MomentumInBaseEvaluable : public levi::DefaultEvaluable {

    ExpressionsServer* m_expressionsServer;
    levi::Expression m_thisExpression;
    levi::Expression m_relativeJacobian;
    levi::Expression m_adjoint;

    levi::Variable m_jointsVariable, m_jointsVelocityVariable, m_baseTwist;

public:

    MomentumInBaseEvaluable(ExpressionsServer* expressionsServer,
                            const levi::Variable &baseTwist)
        : levi::DefaultEvaluable(6,1, "h_B")
          , m_expressionsServer(expressionsServer)
          , m_jointsVariable(expressionsServer->jointsPosition())
          , m_jointsVelocityVariable(expressionsServer->jointsVelocity())
          , m_baseTwist(baseTwist)
    {
        assert(expressionsServer);
        std::string baseFrame = expressionsServer->getFloatingBase();
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::LinkIndex baseIndex = model.getLinkIndex(baseFrame);
        assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

        iDynTree::LinkConstPtr baseLink = model.getLink(baseIndex);
        levi::Constant baseInertia(iDynTree::toEigen(baseLink->getInertia().asMatrix()),"I_b");

        m_thisExpression = baseInertia * baseTwist;
        std::string linkName;

        levi::Expression adjointWrench, adjoint, absoluteVelocity;

        for (iDynTree::LinkIndex l = 0; l < static_cast<int>(model.getNrOfLinks()); ++l) {
            if (l != baseIndex) {
                linkName = model.getLinkName(l);
                m_thisExpression = m_thisExpression + expressionsServer->linkInertiaInBase(l) * (expressionsServer->absoluteVelocity(linkName, baseTwist));
            }
        }

        addDependencies(m_thisExpression);
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        m_evaluationBuffer = m_thisExpression.evaluate();

        return m_evaluationBuffer;
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_thisExpression.clearDerivativesCache();
        m_relativeJacobian.clearDerivativesCache();
        m_adjoint.clearDerivativesCache();
    }

    virtual levi::Expression getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::Expression DynamicalPlanner::Private::MomentumInBaseEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
    assert(column == 0);
    levi::unused(column);

    if (variable->dimension() == m_jointsVariable.rows() && variable->variableName() == m_jointsVariable.name()) {
        return MomentumInBaseExpressionJointsDerivativeExpression(m_expressionsServer, m_baseTwist);
    } else if (m_baseTwist.isDependentFrom(variable)) {
        return MomentumInBaseExpressionBaseTwistDerivativeExpression(m_expressionsServer, m_baseTwist.getColumnDerivative(0, variable));
    } else {
        return m_thisExpression.getColumnDerivative(0, variable);
    }
}

/*-------------------------- Momentum Joints Derivative ---------------------------------------------*/


class DynamicalPlanner::Private::MomentumInBaseJointsDerivativeEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    ExpressionsServer* m_expressionsServer;

    levi::Variable m_baseTwist;

    struct JointsInfo {
        levi::Expression adjointWrenchTimesMotionVector;
        levi::Expression motionVectorTimesChildVelocity;
    };

    std::unordered_map<size_t, JointsInfo> m_jointsInfo;

public:

    MomentumInBaseJointsDerivativeEvaluable(ExpressionsServer* expressionsServer,
                                            const levi::Variable &baseTwist)
        : levi::DefaultEvaluable(6, expressionsServer->jointsPosition().rows(), "d(h_B)/dq")
          , m_expressionsServer(expressionsServer)
          , m_baseTwist(baseTwist)
    {

        assert(expressionsServer);
        std::string baseFrame = expressionsServer->getFloatingBase();
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(expressionsServer->getFloatingBase());
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex baseLink = model.getFrameLink(baseFrameIndex);
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
        assert(ok);

        m_cols.resize(model.getNrOfDOFs(), levi::Null(6,1));

        iDynTree::LinkConstPtr linkPtr;
        iDynTree::IJointConstPtr jointPtr;
        iDynTree::LinkIndex visitedLink, childLink, parentLink;
        size_t jointIndex;

        for (iDynTree::LinkIndex l = 0; l < static_cast<int>(model.getNrOfLinks()); ++l) {

            assert(model.isValidLinkIndex(l));
            linkPtr = model.getLink(l);
            visitedLink = l;

            levi::Expression linkMomentum = expressionsServer->linkInertia(l) * expressionsServer->absoluteVelocity(model.getLinkName(l), baseTwist);

            while(visitedLink != baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                JointsInfo& jointInfo = m_jointsInfo[jointIndex];

                if (!jointInfo.motionVectorTimesChildVelocity.isValidExpression()) {
                    jointInfo.motionVectorTimesChildVelocity = expressionsServer->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), parentLink, childLink) *
                        expressionsServer->absoluteVelocity(model.getLinkName(childLink), baseTwist);

                    jointInfo.adjointWrenchTimesMotionVector = expressionsServer->adjointTransformWrench(baseFrame, model.getLinkName(childLink)) *
                        expressionsServer->motionSubSpaceAsCrossProductWrench(jointPtr->getIndex(), parentLink, childLink);
                }

                levi::Expression adjointWrenchDerivative = jointInfo.adjointWrenchTimesMotionVector *
                    expressionsServer->adjointTransformWrench(model.getLinkName(childLink), model.getLinkName(l)) *
                    linkMomentum;

                levi::Expression velocityDerivative = expressionsServer->linkInertiaInBase(l) *
                    expressionsServer->adjointTransform(model.getLinkName(l), model.getLinkName(childLink)) *
                    jointInfo.motionVectorTimesChildVelocity;

                m_cols[jointIndex] = m_cols[jointIndex] + adjointWrenchDerivative - velocityDerivative;

                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        }

        addDependencies(expressionsServer->jointsPosition(), expressionsServer->jointsVelocity(), baseTwist);

    }

    virtual levi::ColumnExpression col(Eigen::Index col) final {
        return m_cols[static_cast<size_t>(col)];
    }

    virtual levi::ScalarExpression element(Eigen::Index row, Eigen::Index col) final {
        return m_cols[static_cast<size_t>(col)](row, 0);
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        for (size_t j = 0; j < m_cols.size(); ++j) {
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

        for (auto& expression : m_jointsInfo) {
            expression.second.adjointWrenchTimesMotionVector.clearDerivativesCache();
            expression.second.motionVectorTimesChildVelocity.clearDerivativesCache();
        }
    }

};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
MomentumInBaseJointsDerivativeEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->dimension() == m_expressionsServer->jointsPosition().rows() &&
        variable->variableName() == m_expressionsServer->jointsPosition().name()) {
        return m_expressionsServer->momentumInBaseJointsDoubleDerivative(m_baseTwist, column);
    } else {
        return m_cols[static_cast<size_t>(column)].getColumnDerivative(0, variable);
    }
}

/*-------------------------- Momentum Base Twist Derivative ---------------------------------------------*/


class DynamicalPlanner::Private::MomentumInBaseBaseTwistDerivativeEvaluable : public levi::DefaultEvaluable {

    levi::Expression m_thisExpression;

public:

    MomentumInBaseBaseTwistDerivativeEvaluable(ExpressionsServer* expressionsServer,
                                               const levi::Expression &baseTwistJacobian)
        : levi::DefaultEvaluable(6, baseTwistJacobian.cols(), "d(h_B)/d.nu_base * " + baseTwistJacobian.name())
    {

        assert(expressionsServer);

        m_thisExpression = expressionsServer->compositeRigidBodyInertia() * baseTwistJacobian;

        addDependencies(m_thisExpression);

    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        m_evaluationBuffer = m_thisExpression.evaluate();

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_thisExpression.clearDerivativesCache();
    }

};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
MomentumInBaseBaseTwistDerivativeEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    return m_thisExpression.getColumnDerivative(column, variable);
}


/*-------------------------- Momentum Base Joints Double Derivative ---------------------------------------------*/


class DynamicalPlanner::Private::MomentumInBaseJointsDoubleDerivativeEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    std::unordered_set<size_t> m_nonZeros;

public:

    MomentumInBaseJointsDoubleDerivativeEvaluable(ExpressionsServer* es,
                                                  const levi::Variable &baseTwist, long column)
        : levi::DefaultEvaluable(6, es->jointsPosition().rows(), "d(h_B)/dq")
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

        m_cols.resize(model.getNrOfDOFs(), levi::Null(6,1));

        iDynTree::LinkConstPtr linkPtr;
        iDynTree::IJointConstPtr jointPtr;
        iDynTree::LinkIndex visitedLink, childLink, parentLink;
        size_t jointIndex;

        iDynTree::LinkIndex columnChild = traversal.getChildLinkIndexFromJointIndex(model, static_cast<iDynTree::JointIndex>(column));
        std::string columnChildName = model.getLinkName(columnChild);
        iDynTree::LinkIndex columnParent = traversal.getParentLinkIndexFromJointIndex(model, static_cast<iDynTree::JointIndex>(column));

        levi::Expression thisSWrench = es->motionSubSpaceAsCrossProductWrench(static_cast<iDynTree::JointIndex>(column), columnParent, columnChild);

        levi::Expression thisColumnExpressedInC = levi::Null(6,1);

        std::unordered_set<size_t> otherJointsInSameBranch;

        for (iDynTree::LinkIndex l = 0; l < static_cast<int>(model.getNrOfLinks()); ++l) {

            assert(model.isValidLinkIndex(l));
            linkPtr = model.getLink(l);
            visitedLink = l;

            std::vector<size_t> visitedJoints;

            while(visitedLink != baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                if (jointIndex == static_cast<size_t>(column)) {

                    levi::Expression inertiaInC = es->adjointTransformWrench(columnChildName, model.getLinkName(l)) * es->linkInertia(l);

                    levi::Expression motionSubspaceInC = es->adjointTransform(model.getLinkName(l), columnChildName) *
                        es->motionSubSpaceAsCrossProduct(static_cast<iDynTree::JointIndex>(column), columnParent, columnChild) *
                        es->adjointTransform(columnChildName, model.getLinkName(l));

                    levi::Expression linkVel = es->absoluteVelocity(model.getLinkName(l), baseTwist);

                    thisColumnExpressedInC = thisColumnExpressedInC + thisSWrench * inertiaInC * linkVel +
                        inertiaInC * linkVel.getColumnDerivative(0, es->jointsPosition()).col(column);

                    while (visitedLink != baseLink) {
                        jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                        jointIndex = static_cast<size_t>(jointPtr->getIndex());

                        childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
                        parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                        m_cols[jointIndex] = m_cols[jointIndex] +
                            (thisSWrench * inertiaInC - inertiaInC * motionSubspaceInC) * linkVel.getColumnDerivative(0, es->jointsPosition()).col(static_cast<long>(jointIndex));

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

        levi::Expression b_X_c = es->adjointTransformWrench(baseFrame, columnChildName);

        for (size_t modifiedCols : m_nonZeros) {
            jointPtr = model.getJoint(static_cast<iDynTree::JointIndex>(modifiedCols));

            childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
            parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            levi::Expression b_X_c_derivative = es->adjointTransformWrench(baseFrame, model.getLinkName(childLink)) *
                es->motionSubSpaceAsCrossProductWrench(jointPtr->getIndex(), parentLink, childLink) *
                es->adjointTransformWrench(model.getLinkName(childLink), columnChildName);

            m_cols[modifiedCols] = b_X_c_derivative * thisColumnExpressedInC + b_X_c * m_cols[modifiedCols];
        }

        m_nonZeros.insert(otherJointsInSameBranch.begin(), otherJointsInSameBranch.end());

        for (size_t otherCols : otherJointsInSameBranch) {
            m_cols[otherCols] = es->momentumInBaseJointsDoubleDerivative(baseTwist, static_cast<long>(otherCols)).col(column); //here we are exploiting the fact that d_sj*d_sk is equal to
                                                                                                                               //d_sk*d_sj. Thus we compute only the derivatives for which
                                                                                                                               //s_k is closer to the base than s_j
        }

        addDependencies(es->jointsPosition(), es->jointsVelocity(), baseTwist);

    }

    virtual levi::ColumnExpression col(Eigen::Index col) final {
        return m_cols[static_cast<size_t>(col)];
    }

    virtual levi::ScalarExpression element(Eigen::Index row, Eigen::Index col) final {
        return m_cols[static_cast<size_t>(col)](row, 0);
    }

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
MomentumInBaseJointsDoubleDerivativeEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    return m_cols[static_cast<size_t>(column)].getColumnDerivative(0, variable);
}

levi::Expression DynamicalPlanner::Private::MomentumInBaseExpression(ExpressionsServer *expressionsServer, const levi::Variable &baseTwist)
{
    return levi::ExpressionComponent<MomentumInBaseEvaluable>(expressionsServer, baseTwist);
}

levi::Expression DynamicalPlanner::Private::MomentumInBaseExpressionJointsDerivativeExpression(ExpressionsServer *expressionsServer, const levi::Variable &baseTwist)
{
    return levi::ExpressionComponent<MomentumInBaseJointsDerivativeEvaluable>(expressionsServer, baseTwist);
}

levi::Expression DynamicalPlanner::Private::MomentumInBaseExpressionBaseTwistDerivativeExpression(ExpressionsServer *expressionsServer, const levi::Expression& baseTwistJacobian)
{
    return levi::ExpressionComponent<MomentumInBaseBaseTwistDerivativeEvaluable>(expressionsServer, baseTwistJacobian);
}

levi::Expression DynamicalPlanner::Private::MomentumInBaseExpressionJointsDoubleDerivativeExpression(ExpressionsServer *expressionsServer, const levi::Variable &baseTwist, long column)
{
    return levi::ExpressionComponent<MomentumInBaseJointsDoubleDerivativeEvaluable>(expressionsServer, baseTwist, column);
}
