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
    }
}

using namespace DynamicalPlanner::Private;

using AdjointExpression = levi::ExpressionComponent<AdjointTransformEvaluable>;

class DynamicalPlanner::Private::AdjointTransformEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    std::shared_ptr<TimelySharedKinDynComputations> m_timelySharedKinDyn;
    RobotState* m_robotState;
    levi::ScalarVariable m_timeVariable;
    iDynTree::Transform m_transform;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;
    iDynTree::LinkIndex m_baseLink;
    iDynTree::Traversal m_traversal;

    std::vector<levi::Expression> m_jointsDerivative;
    std::vector<levi::ExpressionComponent<
        typename levi::Evaluable<typename levi::DefaultEvaluable::derivative_evaluable::col_type>>> m_colsDerivatives;

public:

    AdjointTransformEvaluable(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                              RobotState *robotState, const std::string &baseFrame, const std::string &targetFrame,
                              levi::Variable jointsVariable, levi::ScalarVariable timeVariable)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (jointsVariable, 6,6, sharedKinDyn->getFloatingBase() + "_X_" + targetFrame)
        , m_timelySharedKinDyn(sharedKinDyn)
        , m_robotState(robotState)
        , m_timeVariable(timeVariable)
    {
        const iDynTree::Model& model = sharedKinDyn->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);
        m_baseLink = model.getFrameLink(m_baseFrame);
        assert(m_baseLink != iDynTree::LINK_INVALID_INDEX);

        bool ok = model.computeFullTreeTraversal(m_traversal, m_baseLink);
        assert(ok);

        m_jointsDerivative.resize(static_cast<size_t>(jointsVariable.rows()), levi::Null(6,6));
        m_colsDerivatives.resize(static_cast<size_t>(jointsVariable.rows()), levi::Null(6,1));

        iDynTree::IJointConstPtr jointPtr;
        size_t jointIndex;
        iDynTree::LinkIndex childLink, parentLink;
        iDynTree::LinkIndex visitedLink = model.getFrameLink(m_targetFrame);
        Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProduct;

        while (visitedLink != m_baseLink) {
            jointPtr = m_traversal.getParentJointFromLinkIndex(visitedLink);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());

            childLink =  m_traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
            parentLink = m_traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            motionSubSpaceAsCrossProduct = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                               childLink,
                                                                                               parentLink).asCrossProductMatrix());

            m_jointsDerivative[jointIndex] =
                - AdjointExpression(sharedKinDyn, robotState, baseFrame, model.getLinkName(childLink),jointsVariable, timeVariable) *
                levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(jointIndex) + "x") *
                AdjointExpression(sharedKinDyn, robotState, model.getLinkName(childLink), targetFrame, jointsVariable, timeVariable);

            visitedLink = m_traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
        }

    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_timelySharedKinDyn->get(m_timeVariable.evaluate());

        m_transform = kinDyn->getRelativeTransform(*m_robotState, m_baseFrame, m_targetFrame);

        m_evaluationBuffer = iDynTree::toEigen(m_transform.asAdjointTransform());

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
AdjointTransformEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name()) { //m_expression contains the jointsVariable specified in the constructor

        for (size_t i = 0; i < m_jointsDerivative.size(); ++i) {
            m_colsDerivatives[i] = m_jointsDerivative[i].col(column);
        }

        return levi::Expression::ComposeByCols(m_colsDerivatives, "d" + name() + "/dq");

    } else {
        return levi::Null(6, variable->dimension());
    }

}


levi::Expression DynamicalPlanner::Private::AdjointTransformExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                                       RobotState *robotState, const std::string &baseFrame,
                                                                       const std::string &targetFrame, levi::Variable jointsVariable,
                                                                       levi::ScalarVariable timeVariable)
{
    return levi::ExpressionComponent<AdjointTransformEvaluable>(sharedKinDyn, robotState, baseFrame,
                                                                targetFrame, jointsVariable, timeVariable);
}


