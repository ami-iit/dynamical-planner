/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/ContactFrictionConstraint.h>
#include <iDynTree/Core/MatrixDynSize.h>

using namespace DynamicalPlanner::Private;

class ContactFrictionConstraint::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    double frictionCoefficient;

    iDynTree::IndexRange forcePointRange;

    iDynTree::Vector3 pointForce;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;
};



ContactFrictionConstraint::ContactFrictionConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                     const std::string &footName, size_t contactIndex)
    : iDynTree::optimalcontrol::Constraint (1, "ContactFriction" + footName + std::to_string(contactIndex))
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->frictionCoefficient = 0.3;

    m_isLowerBounded = false;
    m_isUpperBounded = true;

}

bool ContactFrictionConstraint::setFrictionCoefficient(double frictionCoefficient)
{
    if (frictionCoefficient <= 0.0)
        return false;

    m_pimpl->frictionCoefficient = frictionCoefficient;
    return true;
}

bool ContactFrictionConstraint::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;

    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);

    constraint(0) = m_pimpl->pointForce(0) * m_pimpl->pointForce(0) + m_pimpl->pointForce(1) * m_pimpl->pointForce(1) -
            m_pimpl->frictionCoefficient * m_pimpl->pointForce(2);

    return true;

}

bool ContactFrictionConstraint::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;

    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);

    unsigned int col = static_cast<unsigned int>(m_pimpl->forcePointRange.offset);

    m_pimpl->stateJacobianBuffer(0, col) = 2 * m_pimpl->pointForce(0);
    m_pimpl->stateJacobianBuffer(0, col+1) = 2 * m_pimpl->pointForce(1);
    m_pimpl->stateJacobianBuffer(0, col+2) = -m_pimpl->frictionCoefficient;

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;
}

bool ContactFrictionConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;

    return true;
}

size_t ContactFrictionConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ContactFrictionConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}
