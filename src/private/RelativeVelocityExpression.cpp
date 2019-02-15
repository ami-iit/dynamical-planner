/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeJacobianExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeVelocityExpression.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>


namespace DynamicalPlanner {
    namespace Private {
        class RelativeLeftVelocityEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::RelativeLeftVelocityEvaluable : public levi::DefaultEvaluable {

    std::shared_ptr<TimelySharedKinDynComputations> m_timelySharedKinDyn;
    RobotState* m_robotState;
    std::string m_baseFrameName, m_targetFrameName;
    levi::Variable m_jointsVariable;
    levi::Variable m_jointsVelocityVariable;
    levi::ScalarVariable m_timeVariable;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;
    levi::Expression m_thisExpression;

public:

    RelativeLeftVelocityEvaluable(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                        RobotState *robotState,
                                                        const std::string& baseFrame,
                                                        const std::string &targetFrame,
                                                        levi::Variable jointsVariable,
                                                        levi::Variable jointsVelocityVariable,
                                                        levi::ScalarVariable timeVariable)
        : levi::DefaultEvaluable(6,1, targetFrame + "_V_" + baseFrame + "," + targetFrame)
          , m_timelySharedKinDyn(sharedKinDyn)
          , m_robotState(robotState)
          , m_baseFrameName(baseFrame)
          , m_targetFrameName(targetFrame)
          , m_jointsVariable(jointsVariable)
          , m_jointsVelocityVariable(jointsVelocityVariable)
          , m_timeVariable(timeVariable)
    {
        const iDynTree::Model& model = sharedKinDyn->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);
        iDynTree::LinkIndex baseLink = model.getFrameLink(m_baseFrame);
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);

        levi::Expression jacobian = RelativeLeftJacobianExpression(sharedKinDyn, robotState, baseFrame, targetFrame, jointsVariable, timeVariable);
        m_thisExpression = jacobian * m_jointsVelocityVariable;
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {
        SharedKinDynComputationsPointer kinDyn = m_timelySharedKinDyn->get(m_timeVariable.evaluate());

        iDynTree::Twist velocityInInertial = kinDyn->getFrameVel(*m_robotState, m_targetFrame, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        iDynTree::Twist baseVelocity = kinDyn->getFrameVel(*m_robotState, m_baseFrame, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        iDynTree::Transform relativeTransform = kinDyn->getRelativeTransform(*m_robotState, m_targetFrame, m_baseFrame);

        m_evaluationBuffer = iDynTree::toEigen(velocityInInertial - relativeTransform*baseVelocity);

        return m_evaluationBuffer;
    }

    virtual bool isNew(size_t callerID) final {
        if (m_jointsVariable.isNew() || m_jointsVelocityVariable.isNew()) {
            resetEvaluationRegister();
        }

        return !m_evaluationRegister[callerID];
    }

    virtual levi::ExpressionComponent<derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {
        assert(column == 0);
        levi::unused(column);
        return m_thisExpression.getColumnDerivative(0, variable);
    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable);

};

bool RelativeLeftVelocityEvaluable::isDependentFrom(std::shared_ptr<levi::VariableBase> variable) {
    return ((variable->variableName() == m_jointsVariable.name() && variable->dimension() == m_jointsVariable.rows()) ||
            (variable->variableName() == m_jointsVelocityVariable.name() && variable->dimension() == m_jointsVelocityVariable.rows()));
}

levi::Expression DynamicalPlanner::Private::RelativeLeftVelocityExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                                           RobotState *robotState, const std::string &baseFrame,
                                                                           const std::string &targetFrame, levi::Variable jointsVariable,
                                                                           levi::Variable jointsVelocityVariable, levi::ScalarVariable timeVariable)
{
    return levi::ExpressionComponent<RelativeLeftVelocityEvaluable>(sharedKinDyn, robotState, baseFrame, targetFrame, jointsVariable, jointsVelocityVariable, timeVariable);
}
