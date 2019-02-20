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
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::AdjointTransformEvaluable : public levi::DefaultEvaluable { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    ExpressionsServer* m_expressionsServer;
    iDynTree::Transform m_transform;
    std::string m_baseFrameName, m_targetFrameName, m_parentFrameName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;
    iDynTree::LinkIndex m_baseLink, m_targetLink, m_parentLink;
    iDynTree::IJointConstPtr m_parentJoint = nullptr;
    iDynTree::Transform m_linkToTargetTransform;
    bool m_isConstant = false;

    std::vector<levi::Expression> m_jointsDerivative;
    levi::Expression m_thisMotionSubspace;
    size_t m_parentJointIndex;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_colsDerivatives;

public:

    AdjointTransformEvaluable(ExpressionsServer* expressionsServer,
                              const std::string &baseFrame, const std::string &targetFrame)
        : levi::DefaultEvaluable(6,6, baseFrame + "_X_" + targetFrame)
        , m_expressionsServer(expressionsServer)
        , m_baseFrameName(baseFrame)
        , m_targetFrameName(targetFrame)
    {
        assert(expressionsServer);
        const iDynTree::Model& model = expressionsServer->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);
        m_baseLink = model.getFrameLink(m_baseFrame);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, m_baseLink);
        assert(ok);

        m_jointsDerivative.resize(static_cast<size_t>(expressionsServer->jointsPosition().rows()), levi::Null(6,6));
        m_colsDerivatives.resize(static_cast<size_t>(expressionsServer->jointsPosition().rows()), levi::Null(6,1));

        m_targetLink = model.getFrameLink(m_targetFrame);

        if (m_baseLink != m_targetLink) {

            size_t jointIndex;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = m_targetLink;
            Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProduct;

            m_linkToTargetTransform = model.getFrameTransform(m_targetFrame);

            m_parentJoint = traversal.getParentJointFromLinkIndex(m_targetLink);
            assert(m_parentJoint);
            m_parentJointIndex = static_cast<size_t>(m_parentJoint->getIndex());
            m_parentLink = traversal.getParentLinkIndexFromJointIndex(model, m_parentJoint->getIndex());
            m_parentFrameName = model.getLinkName(m_parentLink);
            motionSubSpaceAsCrossProduct = iDynTree::toEigen(m_parentJoint->getMotionSubspaceVector(0,
                                                                                                    m_targetLink,
                                                                                                    m_parentLink).asCrossProductMatrix());

            m_thisMotionSubspace = levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(m_parentJointIndex) + "x");

            iDynTree::IJointConstPtr jointPtr;
            visitedLink = m_parentLink;
            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                                   visitedLink,
                                                                                                   parentLink).asCrossProductMatrix());

                m_jointsDerivative[jointIndex] = m_expressionsServer->adjointTransform(baseFrame, model.getLinkName(visitedLink)) *
                    levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(jointIndex) + "x") *
                    m_expressionsServer->adjointTransform(model.getLinkName(visitedLink), targetFrame);


                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        } else {
            m_isConstant = true;

            m_evaluationBuffer = iDynTree::toEigen((model.getFrameTransform(m_baseFrame).inverse() * model.getFrameTransform(m_targetFrame)).asAdjointTransform());
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

    virtual bool isNew(size_t callerID) final {
        if (!m_isConstant && (m_expressionsServer->jointsPosition().isNew())) {
            resetEvaluationRegister();
        }

        return !m_evaluationRegister[callerID];
    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) final {
        return (!m_isConstant && (variable->variableName() == m_expressionsServer->jointsPosition().name() && variable->dimension() == m_expressionsServer->jointsPosition().rows()));
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (!m_isConstant && (variable->variableName() == m_expressionsServer->jointsPosition().name() && variable->dimension() == m_expressionsServer->jointsPosition().rows())) {

        m_jointsDerivative[m_parentJointIndex] =
            /**m_expressionsServer->adjointTransform(m_baseFrameName, m_targetFrameName)*/AdjointTransformExpression(m_expressionsServer, m_baseFrameName, m_targetFrameName) *
            m_thisMotionSubspace;

        for (size_t i = 0; i < m_jointsDerivative.size(); ++i) {
            m_colsDerivatives[i] = m_jointsDerivative[i].col(column);
        }

        return levi::Expression::ComposeByCols(m_colsDerivatives, "d" + name() + "/dq");

    } else {
        return levi::Null(6, variable->dimension());
    }

}


class DynamicalPlanner::Private::AdjointTransformWrenchEvaluable : public levi::DefaultEvaluable { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    ExpressionsServer* m_expressionsServer;
    iDynTree::Transform m_transform;
    std::string m_baseFrameName, m_targetFrameName, m_parentFrameName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;
    iDynTree::LinkIndex m_baseLink, m_targetLink, m_parentLink;
    iDynTree::IJointConstPtr m_parentJoint = nullptr;
    iDynTree::Transform m_linkToTargetTransform;
    bool m_isConstant = false;

    std::vector<levi::Expression> m_jointsDerivative;
    levi::Expression m_thisMotionSubspace;
    size_t m_parentJointIndex;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_colsDerivatives;

public:

    AdjointTransformWrenchEvaluable(ExpressionsServer* expressionsServer,
                              const std::string &baseFrame, const std::string &targetFrame)
        : levi::DefaultEvaluable(6,6, baseFrame + "_X*_" + targetFrame)
          , m_expressionsServer(expressionsServer)
          , m_baseFrameName(baseFrame)
          , m_targetFrameName(targetFrame)
    {
        assert(expressionsServer);
        const iDynTree::Model& model = expressionsServer->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);
        m_baseLink = model.getFrameLink(m_baseFrame);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, m_baseLink);
        assert(ok);

        m_jointsDerivative.resize(static_cast<size_t>(expressionsServer->jointsPosition().rows()), levi::Null(6,6));
        m_colsDerivatives.resize(static_cast<size_t>(expressionsServer->jointsPosition().rows()), levi::Null(6,1));

        m_targetLink = model.getFrameLink(m_targetFrame);

        if (m_baseLink != m_targetLink) {

            size_t jointIndex;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = m_targetLink;
            Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProduct;

            m_linkToTargetTransform = model.getFrameTransform(m_targetFrame);

            m_parentJoint = traversal.getParentJointFromLinkIndex(m_targetLink);
            assert(m_parentJoint);
            m_parentJointIndex = static_cast<size_t>(m_parentJoint->getIndex());
            m_parentLink = traversal.getParentLinkIndexFromJointIndex(model, m_parentJoint->getIndex());
            m_parentFrameName = model.getLinkName(m_parentLink);
            motionSubSpaceAsCrossProduct = iDynTree::toEigen(m_parentJoint->getMotionSubspaceVector(0,
                                                                                                    m_targetLink,
                                                                                                    m_parentLink).asCrossProductMatrix());

            m_thisMotionSubspace = levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(m_parentJointIndex) + "x");

            iDynTree::IJointConstPtr jointPtr;
            visitedLink = m_parentLink;
            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                                   visitedLink,
                                                                                                   parentLink).asCrossProductMatrix());

                m_jointsDerivative[jointIndex] = m_expressionsServer->adjointTransformWrench(baseFrame, model.getLinkName(visitedLink)) *
                    levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(jointIndex) + "x") *
                    m_expressionsServer->adjointTransformWrench(model.getLinkName(visitedLink), targetFrame);


                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        } else {
            m_isConstant = true;

            m_evaluationBuffer = iDynTree::toEigen((model.getFrameTransform(m_baseFrame).inverse() * model.getFrameTransform(m_targetFrame)).asAdjointTransformWrench());
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

    virtual bool isNew(size_t callerID) final {
        if (!m_isConstant && (m_expressionsServer->jointsPosition().isNew())) {
            resetEvaluationRegister();
        }

        return !m_evaluationRegister[callerID];
    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) final {
        return (!m_isConstant && (variable->variableName() == m_expressionsServer->jointsPosition().name() && variable->dimension() == m_expressionsServer->jointsPosition().rows()));
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformWrenchEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (!m_isConstant && (variable->variableName() == m_expressionsServer->jointsPosition().name() && variable->dimension() == m_expressionsServer->jointsPosition().rows())) {

        m_jointsDerivative[m_parentJointIndex] =
            /**m_expressionsServer->adjointTransform(m_baseFrameName, m_targetFrameName)*/AdjointTransformWrenchExpression(m_expressionsServer, m_baseFrameName, m_targetFrameName) *
            m_thisMotionSubspace;

        for (size_t i = 0; i < m_jointsDerivative.size(); ++i) {
            m_colsDerivatives[i] = m_jointsDerivative[i].col(column);
        }

        return levi::Expression::ComposeByCols(m_colsDerivatives, "d" + name() + "/dq");

    } else {
        return levi::Null(6, variable->dimension());
    }

}

levi::Expression DynamicalPlanner::Private::AdjointTransformExpression(ExpressionsServer* expressionServer, const std::string &baseFrame,
                                                                       const std::string &targetFrame)
{
    return levi::ExpressionComponent<AdjointTransformEvaluable>(expressionServer, baseFrame, targetFrame);
}

levi::Expression DynamicalPlanner::Private:: AdjointTransformWrenchExpression(ExpressionsServer *expressionServer, const std::string &baseFrame,
                                                                             const std::string &targetFrame)
{
    return levi::ExpressionComponent<AdjointTransformWrenchEvaluable>(expressionServer, baseFrame, targetFrame);
}
