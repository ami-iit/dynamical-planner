/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionErrorExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>


namespace DynamicalPlanner {
    namespace Private {
        class QuaternionErrorEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::QuaternionErrorEvaluable : public levi::DefaultEvaluable {

    std::string m_desiredFrameName;
    ExpressionsServer* m_expressionsServer;
    levi::Variable m_desiredQuaternion;
    iDynTree::FrameIndex m_desiredFrameIndex;
    iDynTree::Rotation m_desiredRotation;
    std::string m_baseFrame;

public:

    QuaternionErrorEvaluable(const std::string &desiredFrame, ExpressionsServer* expressionsServer, const levi::Variable& desiredQuaternion)
        : levi::DefaultEvaluable(4, 1, "q_err_"  +desiredFrame)
          , m_desiredFrameName(desiredFrame)
          , m_expressionsServer(expressionsServer)
          , m_desiredQuaternion(desiredQuaternion)
    {
        m_desiredFrameIndex = expressionsServer->model().getFrameIndex(desiredFrame);
        m_baseFrame = expressionsServer->getFloatingBase();
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        const RobotState& currentState = m_expressionsServer->currentState();
        iDynTree::Transform frameTransform =  m_expressionsServer->currentKinDyn()->getWorldTransform(currentState, m_desiredFrameIndex);

        iDynTree::Vector4 desiredQuaternion_iDyn;
        iDynTree::toEigen(desiredQuaternion_iDyn) = m_desiredQuaternion.evaluate();
        m_desiredRotation.fromQuaternion(desiredQuaternion_iDyn);


        m_evaluationBuffer = iDynTree::toEigen(ErrorQuaternion(frameTransform.getRotation(), m_desiredRotation));


        return m_evaluationBuffer;
    }

    virtual bool isNew(size_t callerID) final {
        if (m_expressionsServer->jointsPosition()->isNew() || m_expressionsServer->baseQuaternion()->isNew() || m_desiredQuaternion.isNew()) {
            resetEvaluationRegister();
        }

        return !m_evaluationRegister[callerID];
    }

    virtual levi::ExpressionComponent<derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {
        assert(column == 0);
        levi::unused(column);

        if (variable->variableName() == m_expressionsServer->jointsPosition()->name() && variable->dimension() == m_expressionsServer->jointsPosition()->rows()) {
            return 0.5 * G_Expression(m_expressionsServer->quaternionError(m_desiredFrameName, m_desiredQuaternion)->asVariable()).transpose() *
                m_expressionsServer->relativeLeftJacobian(m_baseFrame, m_desiredFrameName)->block(3, 0, 3, m_expressionsServer->jointsPosition()->rows());
        }

        if (variable->variableName() == m_expressionsServer->baseQuaternion()->name() && variable->dimension() == m_expressionsServer->baseQuaternion()->rows()) {
            return G_Expression(m_expressionsServer->quaternionError(m_desiredFrameName, m_desiredQuaternion)->asVariable()).transpose() *
                *m_expressionsServer->relativeRotation(m_desiredFrameName, m_baseFrame) * G_Expression(m_expressionsServer->normalizedBaseQuaternion()->asVariable()) *
                m_expressionsServer->normalizedBaseQuaternion()->getColumnDerivative(0, *m_expressionsServer->baseQuaternion());
        }

        if (variable->variableName() == m_desiredQuaternion.name() && variable->dimension() == m_desiredQuaternion.rows()) {
            return levi::ExpressionComponent<derivative_evaluable>(); //this derivative is not provided
        }

        return levi::Null(4, variable->dimension());

    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable);

};

bool QuaternionErrorEvaluable::isDependentFrom(std::shared_ptr<levi::VariableBase> variable) {
    return ((variable->variableName() == m_expressionsServer->jointsPosition()->name() && variable->dimension() == m_expressionsServer->jointsPosition()->rows()) ||
            (variable->variableName() == m_expressionsServer->baseQuaternion()->name() && variable->dimension() == m_expressionsServer->baseQuaternion()->rows()) ||
            (variable->variableName() == m_desiredQuaternion.name() && variable->dimension() == m_desiredQuaternion.rows()));
}

levi::Expression DynamicalPlanner::Private::QuaternionError(const std::string &desiredFrame, ExpressionsServer *expressionsServer, const levi::Variable &desiredQuaternion)
{
    return levi::ExpressionComponent<QuaternionErrorEvaluable>(desiredFrame, expressionsServer, desiredQuaternion);
}
