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

class DynamicalPlanner::Private::AdjointTransformEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    ExpressionsServer* m_expressionsServer;
    iDynTree::Transform m_transform;
    std::string m_baseFrameName, m_targetFrameName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;
    iDynTree::LinkIndex m_baseLink, m_targetLink;

    std::vector<levi::Expression> m_jointsDerivative;
    levi::Expression m_thisMotionSubspace;
    size_t m_parentJointIndex;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_colsDerivatives;

public:

    AdjointTransformEvaluable(ExpressionsServer* expressionsServer,
                              const std::string &baseFrame, const std::string &targetFrame)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (*expressionsServer->jointsPosition(), 6,6, baseFrame + "_X_" + targetFrame)
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

        m_jointsDerivative.resize(static_cast<size_t>(expressionsServer->jointsPosition()->rows()), levi::Null(6,6));
        m_colsDerivatives.resize(static_cast<size_t>(expressionsServer->jointsPosition()->rows()), levi::Null(6,1));

        m_targetLink = model.getFrameLink(m_targetFrame);

        if (m_baseLink != m_targetLink) {

            iDynTree::IJointConstPtr jointPtr;
            size_t jointIndex;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = m_targetLink;
            Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProduct;

            jointPtr = traversal.getParentJointFromLinkIndex(m_targetLink);
            assert(jointPtr);
            m_parentJointIndex = static_cast<size_t>(jointPtr->getIndex());
            parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                               m_targetLink,
                                                                                               parentLink).asCrossProductMatrix());

            m_thisMotionSubspace = levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(m_parentJointIndex) + "x");

            visitedLink = parentLink;
            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                                   visitedLink,
                                                                                                   parentLink).asCrossProductMatrix());

                m_jointsDerivative[jointIndex] = *m_expressionsServer->adjointTransform(baseFrame, model.getLinkName(visitedLink)) *
                    levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(jointIndex) + "x") *
                    *m_expressionsServer->adjointTransform(model.getLinkName(visitedLink), targetFrame);


                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        }

    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_expressionsServer->currentKinDyn();

        m_transform = kinDyn->getRelativeTransform(m_expressionsServer->currentState(), m_baseFrame, m_targetFrame);

        m_evaluationBuffer = iDynTree::toEigen(m_transform.asAdjointTransform());

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name() && (m_baseLink != m_targetLink)) { //m_expression contains the jointsVariable specified in the constructor

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


class DynamicalPlanner::Private::AdjointTransformWrenchEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    ExpressionsServer* m_expressionsServer;
    iDynTree::Transform m_transform;
    std::string m_baseFrameName, m_targetFrameName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;
    iDynTree::LinkIndex m_baseLink, m_targetLink;

    std::vector<levi::Expression> m_jointsDerivative;
    levi::Expression m_thisMotionSubspace;
    size_t m_parentJointIndex;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_colsDerivatives;

public:

    AdjointTransformWrenchEvaluable(ExpressionsServer* expressionServer,
                                    const std::string& baseFrame,
                                    const std::string &targetFrame)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (*expressionServer->jointsPosition(), 6,6, baseFrame + "_X*_" + targetFrame)
          , m_expressionsServer(expressionServer)
          , m_baseFrameName(baseFrame)
          , m_targetFrameName(targetFrame)
    {
        assert(expressionServer);
        const iDynTree::Model& model = expressionServer->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);
        m_baseLink = model.getFrameLink(m_baseFrame);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, m_baseLink);
        assert(ok);

        m_jointsDerivative.resize(static_cast<size_t>(expressionServer->jointsPosition()->rows()), levi::Null(6,6));
        m_colsDerivatives.resize(static_cast<size_t>(expressionServer->jointsPosition()->rows()), levi::Null(6,1));

        m_targetLink = model.getFrameLink(m_targetFrame);

        if (m_baseLink != m_targetLink) {

            iDynTree::IJointConstPtr jointPtr;
            size_t jointIndex;
            iDynTree::LinkIndex parentLink;
            iDynTree::LinkIndex visitedLink = m_targetLink;
            Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProduct;

            jointPtr = traversal.getParentJointFromLinkIndex(m_targetLink);
            assert(jointPtr);
            m_parentJointIndex = static_cast<size_t>(jointPtr->getIndex());
            parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                               m_targetLink,
                                                                                               parentLink).asCrossProductMatrixWrench());

            m_thisMotionSubspace = levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(m_parentJointIndex) + "x");

            visitedLink = parentLink;
            while (visitedLink != m_baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                assert(jointPtr);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                                   visitedLink,
                                                                                                   parentLink).asCrossProductMatrixWrench());

                m_jointsDerivative[jointIndex] = *m_expressionsServer->adjointTransformWrench(baseFrame, model.getLinkName(visitedLink)) *
                    levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(jointIndex) + "x") *
                    *m_expressionsServer->adjointTransformWrench(model.getLinkName(visitedLink), targetFrame);

                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        }

    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_expressionsServer->currentKinDyn();

        m_transform = kinDyn->getRelativeTransform( m_expressionsServer->currentState(), m_baseFrame, m_targetFrame);

        m_evaluationBuffer = iDynTree::toEigen(m_transform.asAdjointTransformWrench());

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformWrenchEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name() && (m_baseLink != m_targetLink)) { //m_expression contains the jointsVariable specified in the constructor

        m_jointsDerivative[m_parentJointIndex] =
            /**m_expressionsServer->adjointTransformWrench(m_baseFrameName, m_targetFrameName)*/ AdjointTransformWrenchExpression(m_expressionsServer, m_baseFrameName, m_targetFrameName) *
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
