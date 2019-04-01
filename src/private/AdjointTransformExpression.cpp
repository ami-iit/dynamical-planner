/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <iDynTree/Model/Model.h>
#include <iDynTree/Model/Traversal.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <cassert>
#include <vector>

namespace DynamicalPlanner {
    namespace Private {
        class AdjointTransformEvaluable;
        class AdjointTransformWrenchEvaluable;
        class AdjointTransformDerivativeEvaluable;
        class AdjointTransformWrenchDerivativeEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::AdjointTransformEvaluable : public levi::DefaultEvaluable {

    ExpressionsServer* m_expressionsServer;
    std::string m_baseFrameName, m_targetFrameName, m_parentFrameName;
    iDynTree::LinkIndex m_targetLink, m_parentLink;
    iDynTree::IJointConstPtr m_parentJoint = nullptr;
    iDynTree::Transform m_linkToTargetTransform;
    bool m_isConstant = false;

    levi::Variable m_jointsVariable;

public:

    AdjointTransformEvaluable(ExpressionsServer* expressionsServer,
                              const std::string &baseFrame, const std::string &targetFrame)
        : levi::DefaultEvaluable(6,6, baseFrame + "_X_" + targetFrame)
        , m_expressionsServer(expressionsServer)
        , m_baseFrameName(baseFrame)
        , m_targetFrameName(targetFrame)
        , m_jointsVariable(expressionsServer->jointsPosition())
    {
        assert(expressionsServer);
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(baseFrame);
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::FrameIndex targetFrameIndex = model.getFrameIndex(targetFrame);
        assert(targetFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex baseLink = model.getFrameLink(baseFrameIndex);
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
        assert(ok);


        m_targetLink = model.getFrameLink(targetFrameIndex);

        if (baseLink != m_targetLink) {

            m_linkToTargetTransform = model.getFrameTransform(targetFrameIndex);
            m_parentJoint = traversal.getParentJointFromLinkIndex(m_targetLink);
            assert(m_parentJoint);
            m_parentLink = traversal.getParentLinkIndexFromJointIndex(model, m_parentJoint->getIndex());
            m_parentFrameName = model.getLinkName(m_parentLink);

        } else {
            m_isConstant = true;

            m_evaluationBuffer = iDynTree::toEigen((model.getFrameTransform(baseFrameIndex).inverse() * model.getFrameTransform(targetFrameIndex)).asAdjointTransform());
        }

        if (!m_isConstant) {
            addDependencies(m_jointsVariable);
        }

    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        if (!m_isConstant) {
            m_evaluationBuffer = m_expressionsServer->adjointTransform(m_baseFrameName, m_parentFrameName).evaluate() *
                iDynTree::toEigen((m_parentJoint->getTransform(m_expressionsServer->currentState().s, m_parentLink, m_targetLink) *
                                   m_linkToTargetTransform).asAdjointTransform());
        }

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformEvaluable::getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (!m_isConstant && (variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows())) {

        return m_expressionsServer->adjointTransformJointsDerivative(m_baseFrameName, m_targetFrameName, column);

    } else {
        return levi::Null(6, variable->dimension());
    }

}


class DynamicalPlanner::Private::AdjointTransformWrenchEvaluable : public levi::DefaultEvaluable {

    ExpressionsServer* m_expressionsServer;
    std::string m_baseFrameName, m_targetFrameName, m_parentFrameName;
    iDynTree::LinkIndex m_targetLink, m_parentLink;
    iDynTree::IJointConstPtr m_parentJoint = nullptr;
    iDynTree::Transform m_linkToTargetTransform;
    bool m_isConstant = false;

    levi::Variable m_jointsVariable;

public:

    AdjointTransformWrenchEvaluable(ExpressionsServer* expressionsServer,
                                    const std::string &baseFrame, const std::string &targetFrame)
        : levi::DefaultEvaluable(6,6, baseFrame + "_X*_" + targetFrame)
          , m_expressionsServer(expressionsServer)
          , m_baseFrameName(baseFrame)
          , m_targetFrameName(targetFrame)
          , m_jointsVariable(expressionsServer->jointsPosition())
    {
        assert(expressionsServer);
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(baseFrame);
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::FrameIndex targetFrameIndex = model.getFrameIndex(targetFrame);
        assert(targetFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex baseLink = model.getFrameLink(baseFrameIndex);
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
        assert(ok);

        m_targetLink = model.getFrameLink(targetFrameIndex);

        if (baseLink != m_targetLink) {

            m_linkToTargetTransform = model.getFrameTransform(targetFrameIndex);
            m_parentJoint = traversal.getParentJointFromLinkIndex(m_targetLink);
            assert(m_parentJoint);
            m_parentLink = traversal.getParentLinkIndexFromJointIndex(model, m_parentJoint->getIndex());
            m_parentFrameName = model.getLinkName(m_parentLink);

        } else {
            m_isConstant = true;

            m_evaluationBuffer = iDynTree::toEigen((model.getFrameTransform(baseFrameIndex).inverse() * model.getFrameTransform(targetFrameIndex)).asAdjointTransformWrench());
        }

        if (!m_isConstant) {
            addDependencies(m_jointsVariable);
        }
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        if (!m_isConstant) {
            m_evaluationBuffer = m_expressionsServer->adjointTransformWrench(m_baseFrameName, m_parentFrameName).evaluate() *
                iDynTree::toEigen((m_parentJoint->getTransform(m_expressionsServer->currentState().s, m_parentLink, m_targetLink) *
                                   m_linkToTargetTransform).asAdjointTransformWrench());
        }

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformWrenchEvaluable::getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (!m_isConstant && (variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows())) {

        return m_expressionsServer->adjointTransformWrenchJointsDerivative(m_baseFrameName, m_targetFrameName, column);

    } else {
        return levi::Null(6, variable->dimension());
    }

}

class DynamicalPlanner::Private::AdjointTransformDerivativeEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::Expression> m_jointsDerivative;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    std::vector<size_t> m_nonZeros;

    levi::Variable m_jointsVariable;

public:

    AdjointTransformDerivativeEvaluable(ExpressionsServer *expressionsServer,
                                        const std::string &baseFrame,
                                        const std::string &targetFrame,
                                        long column)
        : levi::DefaultEvaluable(6, expressionsServer->jointsPosition().rows(), "d(" + baseFrame + "_X_" + targetFrame + ")[:," + std::to_string(column) + "]/dq")
          , m_jointsVariable(expressionsServer->jointsPosition())
    {

        assert(expressionsServer);
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(baseFrame);
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::FrameIndex targetFrameIndex = model.getFrameIndex(targetFrame);
        assert(targetFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex m_baseLink = model.getFrameLink(baseFrameIndex);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, m_baseLink);
        assert(ok);

        m_jointsDerivative.resize(model.getNrOfDOFs(), levi::Null(6,6));
        m_cols.resize(model.getNrOfDOFs(), levi::Null(6,1));

        iDynTree::LinkIndex targetLink = model.getFrameLink(targetFrameIndex);

        if (m_baseLink != targetLink) {

            size_t jointIndex;
            iDynTree::IJointConstPtr jointPtr;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = targetLink;

            jointPtr = traversal.getParentJointFromLinkIndex(targetLink);
            assert(jointPtr);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());
            parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            m_jointsDerivative[jointIndex] =
                expressionsServer->adjointTransform(baseFrame, targetFrame) * expressionsServer->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), parentLink, targetLink);

            m_nonZeros.push_back(jointIndex);

            visitedLink = parentLink;

            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                m_jointsDerivative[jointIndex] = expressionsServer->adjointTransform(baseFrame, model.getLinkName(visitedLink)) *
                    expressionsServer->motionSubSpaceAsCrossProduct(jointPtr->getIndex(), parentLink, visitedLink) *
                    expressionsServer->adjointTransform(model.getLinkName(visitedLink), targetFrame);

                m_nonZeros.push_back(jointIndex);

                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }

            for (size_t i = 0; i < m_jointsDerivative.size(); ++i) {
                m_cols[i] = m_jointsDerivative[i].col(static_cast<Eigen::Index>(column));
            }

        }

        addDependencies(m_jointsVariable);

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

};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformDerivativeEvaluable::getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (m_nonZeros.size() && (variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows())) {

        return m_cols[static_cast<size_t>(column)].getColumnDerivative(0, variable);

    } else {
        return levi::Null(6, variable->dimension());
    }

}

class DynamicalPlanner::Private::AdjointTransformWrenchDerivativeEvaluable : public levi::DefaultEvaluable {

    std::vector<levi::Expression> m_jointsDerivative;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_cols;

    std::vector<size_t> m_nonZeros;

    levi::Variable m_jointsVariable;

public:

    AdjointTransformWrenchDerivativeEvaluable(ExpressionsServer *expressionsServer,
                                              const std::string &baseFrame,
                                              const std::string &targetFrame,
                                              long column)
        : levi::DefaultEvaluable(6, expressionsServer->jointsPosition().rows(), "d(" + baseFrame + "_X*_" + targetFrame + ")[:," + std::to_string(column) + "]/dq")
          , m_jointsVariable(expressionsServer->jointsPosition())
    {

        assert(expressionsServer);
        const iDynTree::Model& model = expressionsServer->model();
        iDynTree::FrameIndex baseFrameIndex = model.getFrameIndex(baseFrame);
        assert(baseFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::FrameIndex targetFrameIndex = model.getFrameIndex(targetFrame);
        assert(targetFrameIndex != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex m_baseLink = model.getFrameLink(baseFrameIndex);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, m_baseLink);
        assert(ok);

        m_jointsDerivative.resize(model.getNrOfDOFs(), levi::Null(6,6));
        m_cols.resize(model.getNrOfDOFs(), levi::Null(6,1));

        iDynTree::LinkIndex targetLink = model.getFrameLink(targetFrameIndex);

        if (m_baseLink != targetLink) {

            size_t jointIndex;
            iDynTree::IJointConstPtr jointPtr;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = targetLink;

            jointPtr = traversal.getParentJointFromLinkIndex(targetLink);
            assert(jointPtr);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());
            parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            m_jointsDerivative[jointIndex] =
                expressionsServer->adjointTransformWrench(baseFrame, targetFrame) * expressionsServer->motionSubSpaceAsCrossProductWrench(jointPtr->getIndex(), parentLink, targetLink);

            m_nonZeros.push_back(jointIndex);

            visitedLink = parentLink;

            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                m_jointsDerivative[jointIndex] = expressionsServer->adjointTransformWrench(baseFrame, model.getLinkName(visitedLink)) *
                    expressionsServer->motionSubSpaceAsCrossProductWrench(jointPtr->getIndex(), parentLink, visitedLink) *
                    expressionsServer->adjointTransformWrench(model.getLinkName(visitedLink), targetFrame);

                m_nonZeros.push_back(jointIndex);

                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }

            for (size_t i = 0; i < m_jointsDerivative.size(); ++i) {
                m_cols[i] = m_jointsDerivative[i].col(static_cast<Eigen::Index>(column));
            }

        }

        addDependencies(m_jointsVariable);

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

};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformWrenchDerivativeEvaluable::getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (m_nonZeros.size() && (variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows())) {

        return m_cols[static_cast<size_t>(column)].getColumnDerivative(0, variable);

    } else {
        return levi::Null(6, variable->dimension());
    }

}

levi::Expression DynamicalPlanner::Private::AdjointTransformExpression(ExpressionsServer* expressionServer, const std::string &baseFrame,
                                                                       const std::string &targetFrame)
{
    return levi::ExpressionComponent<AdjointTransformEvaluable>(expressionServer, baseFrame, targetFrame);
}

levi::Expression DynamicalPlanner::Private::AdjointTransformExpressionJointsDerivative(ExpressionsServer *expressionsServer,
                                                                                       const std::string &baseFrame,
                                                                                       const std::string &targetFrame,
                                                                                       long column)
{
    return levi::ExpressionComponent<AdjointTransformDerivativeEvaluable>(expressionsServer, baseFrame, targetFrame, column);
}

levi::Expression DynamicalPlanner::Private:: AdjointTransformWrenchExpression(ExpressionsServer *expressionServer, const std::string &baseFrame,
                                                                             const std::string &targetFrame)
{
    return levi::ExpressionComponent<AdjointTransformWrenchEvaluable>(expressionServer, baseFrame, targetFrame);
}

levi::Expression DynamicalPlanner::Private::AdjointTransformWrenchExpressionJointsDerivative(ExpressionsServer *expressionsServer,
                                                                                             const std::string &baseFrame,
                                                                                             const std::string &targetFrame,
                                                                                             long column)
{
    return levi::ExpressionComponent<AdjointTransformWrenchDerivativeEvaluable>(expressionsServer, baseFrame, targetFrame, column);
}
