/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AbsoluteVelocityExpression.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>


namespace DynamicalPlanner {
    namespace Private {
        class AbsoluteLeftVelocityEvaluable;
        class AbsoluteLeftVelocityJointsDerivativeEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::AbsoluteLeftVelocityEvaluable : public levi::DefaultEvaluable {

    ExpressionsServer* m_expressionsServer;
    std::string m_targetFrameName;
    iDynTree::FrameIndex m_targetFrame;
    levi::Expression m_thisExpression;
    levi::Expression m_relativeJacobian;
    levi::Expression m_adjoint;

    levi::Variable m_jointsVariable, m_jointsVelocityVariable, m_baseTwist;

public:

    AbsoluteLeftVelocityEvaluable(ExpressionsServer* expressionsServer,
                                  const levi::Variable &baseTwist,
                                  const std::string &targetFrame)
        : levi::DefaultEvaluable(6,1, targetFrame + "_V_A," + targetFrame)
          , m_expressionsServer(expressionsServer)
          , m_targetFrameName(targetFrame)
          , m_jointsVariable(expressionsServer->jointsPosition())
          , m_jointsVelocityVariable(expressionsServer->jointsVelocity())
          , m_baseTwist(baseTwist)
    {
        const iDynTree::Model& model = expressionsServer->model();
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);

        m_thisExpression = expressionsServer->relativeVelocity(expressionsServer->getFloatingBase(), targetFrame) +
            expressionsServer->adjointTransform(targetFrame, expressionsServer->getFloatingBase()) * m_baseTwist;

        m_relativeJacobian = expressionsServer->relativeLeftJacobian(expressionsServer->getFloatingBase(), targetFrame);

        m_adjoint = expressionsServer->adjointTransform(targetFrame, expressionsServer->getFloatingBase());

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

levi::Expression DynamicalPlanner::Private::AbsoluteLeftVelocityEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
    assert(column == 0);
    levi::unused(column);

    if (variable->dimension() == m_jointsVariable.rows() && variable->variableName() == m_jointsVariable.name()) {
        return m_expressionsServer->absoluteVelocityJointsDerivative(m_targetFrameName, m_baseTwist);
    } else if (variable->dimension() == m_jointsVelocityVariable.rows() && variable->variableName() == m_jointsVelocityVariable.name()) {
        return m_relativeJacobian;
    } else if (m_baseTwist.isDependentFrom(variable)) {
        return m_adjoint * m_baseTwist.getColumnDerivative(0, variable);
    } else {
        return levi::Null(6, variable->dimension());
    }

}

class DynamicalPlanner::Private::AbsoluteLeftVelocityJointsDerivativeEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    std::vector<size_t> m_nonZeros;

    levi::Variable m_jointsVariable;

public:

    AbsoluteLeftVelocityJointsDerivativeEvaluable(ExpressionsServer* expressionsServer,
                                                  const levi::Variable &baseTwist,
                                                  const std::string &targetFrame)
        : levi::DefaultEvaluable(6, expressionsServer->jointsPosition().rows(), "d(" + targetFrame + "_V_A," + targetFrame + ")/dq")
          , m_jointsVariable(expressionsServer->jointsPosition())
    {

        assert(expressionsServer);
        std::string baseFrame = expressionsServer->getFloatingBase();
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(expressionsServer->getFloatingBase());
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::FrameIndex targetFrameIndex = model.getFrameIndex(targetFrame);
        assert(targetFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex m_baseLink = model.getFrameLink(baseFrameIndex);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, m_baseLink);
        assert(ok);

        m_cols.resize(model.getNrOfDOFs(), levi::Null(6,1));

        iDynTree::LinkIndex targetLink = model.getFrameLink(targetFrameIndex);

        if (m_baseLink != targetLink) {

            size_t jointIndex;
            iDynTree::IJointConstPtr jointPtr;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = targetLink;

            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                m_cols[jointIndex] = expressionsServer->adjointTransform(targetFrame, model.getLinkName(visitedLink)) *
                    (-expressionsServer->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), parentLink, visitedLink)) *
                    (expressionsServer->absoluteVelocity(model.getLinkName(visitedLink), baseTwist));

                m_nonZeros.push_back(jointIndex);

                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        }

        addDependencies(m_jointsVariable, expressionsServer->jointsVelocity(), baseTwist);

    }

    virtual levi::ColumnExpression col(Eigen::Index col) final {
        return m_cols[static_cast<size_t>(col)];
    }

    virtual levi::ScalarExpression element(Eigen::Index row, Eigen::Index col) final {
        return m_cols[static_cast<size_t>(col)](row, 0);
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        for (size_t nonZero : m_nonZeros) {
            m_evaluationBuffer.col(static_cast<Eigen::Index>(nonZero)) = m_cols[nonZero].evaluate();
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
    }

};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AbsoluteLeftVelocityJointsDerivativeEvaluable::getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    return m_cols[static_cast<size_t>(column)].getColumnDerivative(0, variable);
}

levi::Expression DynamicalPlanner::Private::AbsoluteLeftVelocityExpression(ExpressionsServer *expressionsServer, const levi::Variable &baseTwist,
                                                                           const std::string &targetFrame)
{
    return levi::ExpressionComponent<AbsoluteLeftVelocityEvaluable>(expressionsServer, baseTwist, targetFrame);
}

levi::Expression DynamicalPlanner::Private::AbsoluteLeftVelocityJointsDerivativeExpression(ExpressionsServer *expressionsServer, const levi::Variable &baseTwist,
                                                                                           const std::string &targetFrame)
{
    return levi::ExpressionComponent<AbsoluteLeftVelocityJointsDerivativeEvaluable>(expressionsServer, baseTwist, targetFrame);
}
