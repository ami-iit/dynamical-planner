/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/MomentumInBaseExpression.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>


namespace DynamicalPlanner {
    namespace Private {
        class MomentumInBaseEvaluable;
        class MomentumInBaseJointsDerivativeEvaluable;
        class MomentumInBaseJointsDoubleDerivativeEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

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
                adjoint = expressionsServer->adjointTransform(linkName, baseFrame);
                adjointWrench = expressionsServer->adjointTransformWrench(baseFrame, linkName);
                levi::Constant inertia(iDynTree::toEigen(model.getLink(l)->getInertia().asMatrix()),"I_" + linkName);

                m_thisExpression = m_thisExpression +  (adjointWrench * inertia) * (expressionsServer->absoluteVelocity(linkName, baseTwist.asVariable()));
            }
        }

        addDependencies(m_jointsVariable, m_jointsVelocityVariable, m_baseTwist);
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
    } else {
        return m_thisExpression.getColumnDerivative(0, variable);
    }
}

class DynamicalPlanner::Private::MomentumInBaseJointsDerivativeEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    struct JointsInfo {
        levi::Expression adjointWrenchTimesMotionVector;
        levi::Expression motionVectorTimesChildVelocity;
    };

    struct LinkInfo {

    };

    std::unordered_map<size_t, JointsInfo> m_jointsInfo;

public:

    MomentumInBaseJointsDerivativeEvaluable(ExpressionsServer* expressionsServer,
                                            const levi::Variable &baseTwist)
        : levi::DefaultEvaluable(6, expressionsServer->jointsPosition().rows(), "d(h_B)/dq")
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
            levi::Expression linkInertiaInBase = expressionsServer->adjointTransformWrench(baseFrame, model.getLinkName(l)) * expressionsServer->linkInertia(l);

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

                levi::Expression velocityDerivative = linkInertiaInBase *
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
    getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;

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
MomentumInBaseJointsDerivativeEvaluable::getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
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
