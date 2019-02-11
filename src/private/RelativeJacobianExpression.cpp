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
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeJacobianExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <cassert>
#include <vector>

namespace DynamicalPlanner {
    namespace Private {
        class RelativeJacobianEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::RelativeJacobianEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    std::shared_ptr<TimelySharedKinDynComputations> m_timelySharedKinDyn;
    RobotState* m_robotState;
    levi::ScalarVariable m_timeVariable;
    iDynTree::MatrixDynSize m_jacobian;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;

    std::vector<levi::Expression> m_columns;

public:

    RelativeJacobianEvaluable(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                              RobotState *robotState, const std::string &baseFrame, const std::string &targetFrame,
                              levi::Variable jointsVariable, levi::ScalarVariable timeVariable)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (jointsVariable, 6, jointsVariable.rows(), sharedKinDyn->getFloatingBase() + "_X_" + targetFrame)
          , m_timelySharedKinDyn(sharedKinDyn)
          , m_robotState(robotState)
          , m_timeVariable(timeVariable)
    {
        const iDynTree::Model& model = sharedKinDyn->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex baseLink = model.getFrameLink(m_baseFrame);
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);

        iDynTree::Traversal traversal;

        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
        assert(ok);

        m_columns.resize(static_cast<size_t>(jointsVariable.rows()), levi::Null(6,1));

        iDynTree::IJointConstPtr jointPtr;
        size_t jointIndex;
        iDynTree::LinkIndex childLink, parentLink;
        iDynTree::LinkIndex visitedLink = model.getFrameLink(m_targetFrame);
        Eigen::Matrix<double, 6,1> motionSubSpaceVector;

        while (visitedLink != baseLink) {
            jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());

            childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
            parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

            motionSubSpaceVector = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                       childLink,
                                                                                       parentLink));

            m_columns[jointIndex] = AdjointTransformExpression(sharedKinDyn, robotState, targetFrame, model.getLinkName(childLink), jointsVariable, timeVariable) *
                levi::Constant(motionSubSpaceVector, "s_" + std::to_string(jointIndex));

            visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
        }

    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_timelySharedKinDyn->get(m_timeVariable.evaluate());

        bool ok = kinDyn->getRelativeJacobian(*m_robotState, m_baseFrame, m_targetFrame, m_jacobian, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        assert(ok);

        m_evaluationBuffer = iDynTree::toEigen(m_jacobian);

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
RelativeJacobianEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name()) { //m_expression contains the jointsVariable specified in the constructor

        return m_columns[static_cast<size_t>(column)].getColumnDerivative(0, variable);

    } else {
        return levi::Null(6, variable->dimension());
    }

}


levi::Expression DynamicalPlanner::Private::RelativeJacobianExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                                       RobotState *robotState, const std::string &baseFrame,
                                                                       const std::string &targetFrame, levi::Variable jointsVariable,
                                                                       levi::ScalarVariable timeVariable)
{
    return levi::ExpressionComponent<RelativeJacobianEvaluable>(sharedKinDyn, robotState, baseFrame,
                                                                targetFrame, jointsVariable, timeVariable);
}


