/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/PhantomForcesCost.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::PhantomForcesCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    HyperbolicSecant activation;

    iDynTree::IndexRange positionPointRange, forcePointRange;
    iDynTree::Vector3 pointPosition, pointForce;

    iDynTree::VectorDynSize stateJacobianBuffer, controlJacobianBuffer;
};



PhantomForcesCost::PhantomForcesCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                     const std::string &footName, size_t contactIndex, const HyperbolicSecant &forceActivation)
    : iDynTree::optimalcontrol::Cost ("PhantomForces" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->activation = forceActivation;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();
}

PhantomForcesCost::~PhantomForcesCost()
{ }

bool PhantomForcesCost::costEvaluation(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double fz = m_pimpl->pointForce(2);

    double phantomForce = (1- delta) * fz;

    costValue = 0.5 * phantomForce * phantomForce;

    return true;
}

bool PhantomForcesCost::costFirstPartialDerivativeWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);

    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double deltaDerivative = m_pimpl->activation.evalDerivative(m_pimpl->pointPosition(2));
    double fz = m_pimpl->pointForce(2);
    unsigned int pzCol = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);
    unsigned int fzCol = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    double phantomForce = (1- delta) * fz;

    m_pimpl->stateJacobianBuffer(pzCol) = - phantomForce * deltaDerivative * fz;

    m_pimpl->stateJacobianBuffer(fzCol) = phantomForce * (1- delta);

    partialDerivative = m_pimpl->stateJacobianBuffer;

    return true;
}

bool PhantomForcesCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->controlJacobianBuffer;
    return true;
}

bool PhantomForcesCost::costSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);

    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double deltaDerivative = m_pimpl->activation.evalDerivative(m_pimpl->pointPosition(2));
    double deltaDoubleDerivative = m_pimpl->activation.evalDoubleDerivative(m_pimpl->pointPosition(2));
    double fz = m_pimpl->pointForce(2);
    unsigned int pzCol = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);
    unsigned int fzCol = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    double phantomForce = (1- delta) * fz;

    partialDerivative(pzCol, pzCol) = (deltaDerivative * fz) * (deltaDerivative * fz) - phantomForce * deltaDoubleDerivative;
    partialDerivative(fzCol, fzCol) = (1- delta) * (1- delta);
    partialDerivative(fzCol, pzCol) = -2.0 * phantomForce * deltaDerivative;
    partialDerivative(pzCol, fzCol) = partialDerivative(fzCol, pzCol);

    return true;
}

bool PhantomForcesCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool PhantomForcesCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}
