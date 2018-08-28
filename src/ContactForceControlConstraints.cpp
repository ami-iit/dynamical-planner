/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/ContactForceControlConstraints.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/MatrixFixSize.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class ContactForceControlConstraints::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    HyperbolicSecant activation;
    iDynTree::Vector3 maximumDerivatives;
    iDynTree::Matrix3x3 dissipationRatios;

    iDynTree::IndexRange positionPointRange, forcePointRange, forceControlRange;
    iDynTree::Vector3 pointPosition, pointForce, pointForceControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;
};




ContactForceControlConstraints::ContactForceControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                               const std::string &footName, size_t contactIndex, const HyperbolicSecant &forceActivation,
                                                               const iDynTree::Vector3 &maximumDerivatives, const iDynTree::Vector3 &dissipationRatios)
    : iDynTree::optimalcontrol::Constraint (6, "ForceControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->activation = forceActivation;
    m_pimpl->maximumDerivatives = maximumDerivatives;
    iDynTree::toEigen(m_pimpl->dissipationRatios) = iDynTree::toEigen(dissipationRatios).asDiagonal();

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->forceControlRange = controlVariables.getIndexRange(footName + "ForceControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(6, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(6, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(6);

    m_isLowerBounded = true;
    m_isUpperBounded = false;
}

bool ContactForceControlConstraints::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);
    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));

    iDynTree::toEigen(m_pimpl->constraintValues).topRows<3>() = delta * iDynTree::toEigen(m_pimpl->maximumDerivatives) -
            (1- delta) * iDynTree::toEigen(m_pimpl->dissipationRatios) * iDynTree::toEigen(m_pimpl->pointForce) - iDynTree::toEigen(m_pimpl->pointForceControl);

    iDynTree::toEigen(m_pimpl->constraintValues).bottomRows<3>() = iDynTree::toEigen(m_pimpl->pointForceControl) + delta * iDynTree::toEigen(m_pimpl->maximumDerivatives) +
            (1- delta) * iDynTree::toEigen(m_pimpl->dissipationRatios) * iDynTree::toEigen(m_pimpl->pointForce);

    constraint = m_pimpl->constraintValues;

    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);

    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double deltaDerivative = m_pimpl->activation.evalDerivative(m_pimpl->pointPosition(2));

    iDynTree::toEigen(m_pimpl->stateJacobianBuffer).block<3,1>(0, m_pimpl->positionPointRange.offset + 2) = deltaDerivative * iDynTree::toEigen(m_pimpl->maximumDerivatives) +
            deltaDerivative * iDynTree::toEigen(m_pimpl->dissipationRatios) * iDynTree::toEigen(m_pimpl->pointForce);

    iDynTree::toEigen(m_pimpl->stateJacobianBuffer).block<3,3>(0, m_pimpl->forcePointRange.offset) = -(1- delta) * iDynTree::toEigen(m_pimpl->dissipationRatios);

    iDynTree::toEigen(m_pimpl->stateJacobianBuffer).block<3,1>(3, m_pimpl->positionPointRange.offset + 2) = deltaDerivative * iDynTree::toEigen(m_pimpl->maximumDerivatives) -
            deltaDerivative * iDynTree::toEigen(m_pimpl->dissipationRatios) * iDynTree::toEigen(m_pimpl->pointForce);

    iDynTree::toEigen(m_pimpl->stateJacobianBuffer).block<3,3>(3, m_pimpl->forcePointRange.offset) = (1- delta) * iDynTree::toEigen(m_pimpl->dissipationRatios);

    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<3,3>(0, m_pimpl->forceControlRange.offset).setIdentity();
    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<3,3>(0, m_pimpl->forceControlRange.offset) *= -1;

    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<3,3>(3, m_pimpl->forceControlRange.offset).setIdentity();

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ContactForceControlConstraints::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ContactForceControlConstraints::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}
