/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/QuaternionNormConstraint.h>
#include <DynamicalPlannerPrivate/QuaternionUtils.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class QuaternionNormConstraint::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;
    iDynTree::IndexRange baseQuaternionRange;
    iDynTree::Vector4 baseQuaternion;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;
};


QuaternionNormConstraint::QuaternionNormConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables)
    : iDynTree::optimalcontrol::Constraint (1, "QuaternionNorm")
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_isLowerBounded = true;
    m_isUpperBounded = true;
    m_upperBound.zero();
    m_lowerBound.zero();

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());
    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(m_pimpl->stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();
    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();
}

QuaternionNormConstraint::~QuaternionNormConstraint()
{ }

void QuaternionNormConstraint::setEqualityTolerance(double tolerance)
{
    assert(tolerance > 0);

    iDynTree::toEigen(m_lowerBound).setConstant(-tolerance/2.0);
    iDynTree::toEigen(m_upperBound).setConstant(tolerance/2.0);
}

bool QuaternionNormConstraint::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->baseQuaternion = m_pimpl->stateVariables(m_pimpl->baseQuaternionRange);
    constraint(0) = QuaternionSquaredNorm(m_pimpl->baseQuaternion) - 1.0;
    return true;
}

bool QuaternionNormConstraint::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->baseQuaternion = m_pimpl->stateVariables(m_pimpl->baseQuaternionRange);
    iDynTree::toEigen(m_pimpl->stateJacobianBuffer).block<1, 4>(0, m_pimpl->baseQuaternionRange.offset) = 2 * iDynTree::toEigen(m_pimpl->baseQuaternion).transpose();
    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool QuaternionNormConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t QuaternionNormConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t QuaternionNormConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();

}
