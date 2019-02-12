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
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeQuaternionExpression.h>
#include <cassert>
#include <vector>

namespace DynamicalPlanner {
    namespace Private {
        class RelativeQuaternionEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::RelativeQuaternionEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    std::shared_ptr<TimelySharedKinDynComputations> m_timelySharedKinDyn;
    RobotState* m_robotState;
    levi::ScalarVariable m_timeVariable;
    std::string m_baseName, m_targetName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;

    levi::Expression m_relativeJacobian;

public:

    RelativeQuaternionEvaluable(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                RobotState *robotState, const std::string &baseFrame, const std::string &targetFrame,
                                levi::Variable jointsVariable, levi::ScalarVariable timeVariable, const levi::Expression& relativeJacobian)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (jointsVariable, 4, 1, baseFrame + "_rho_" + targetFrame)
          , m_timelySharedKinDyn(sharedKinDyn)
          , m_robotState(robotState)
          , m_timeVariable(timeVariable)
          , m_baseName(baseFrame)
          , m_targetName(targetFrame)
    {
        const iDynTree::Model& model = sharedKinDyn->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);

        m_relativeJacobian = relativeJacobian.block(3, 0, 3, jointsVariable.rows());
    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_timelySharedKinDyn->get(m_timeVariable.evaluate());

        m_evaluationBuffer = iDynTree::toEigen(kinDyn->getRelativeTransform(*m_robotState, m_baseFrame, m_targetFrame).getRotation().asQuaternion());

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
RelativeQuaternionEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name()) { //m_expression contains the jointsVariable specified in the constructor

        assert(column == 0);

        levi::unused(column);

        levi::Expression thisQuaternion = RelativeQuaternionExpression(m_timelySharedKinDyn, m_robotState, m_baseName, m_targetName, m_expression, m_timeVariable);

        levi::Expression leftQuaternionMap = 0.5 * G_Expression(thisQuaternion.asVariable()).transpose();

        return leftQuaternionMap * m_relativeJacobian;
    } else {
        return levi::Null(4, variable->dimension());
    }

}


levi::Expression DynamicalPlanner::Private::RelativeQuaternionExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                                         RobotState *robotState, const std::string &baseFrame,
                                                                         const std::string &targetFrame, levi::Variable jointsVariable,
                                                                         levi::ScalarVariable timeVariable)
{

    levi::Expression jacobian = RelativeLeftJacobianExpression(sharedKinDyn, robotState, baseFrame, targetFrame, jointsVariable, timeVariable);

    return levi::ExpressionComponent<RelativeQuaternionEvaluable>(sharedKinDyn, robotState, baseFrame,
                                                                  targetFrame, jointsVariable, timeVariable, jacobian);
}

levi::Expression DynamicalPlanner::Private::RelativeQuaternionExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                                         RobotState *robotState, const std::string &baseFrame,
                                                                         const std::string &targetFrame, levi::Variable jointsVariable,
                                                                         levi::ScalarVariable timeVariable, const levi::Expression& relativeJacobian)
{

    return levi::ExpressionComponent<RelativeQuaternionEvaluable>(sharedKinDyn, robotState, baseFrame,
                                                                  targetFrame, jointsVariable, timeVariable, relativeJacobian);
}


