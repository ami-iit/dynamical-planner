/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Constraints/DynamicalComplementarityConstraint.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>

using namespace DynamicalPlanner::Private;

class DynamicalComplementarityConstraint::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    double dissipationGain;

    iDynTree::IndexRange positionPointRange, velocityControlRange, forcePointRange, forceControlRange;
    iDynTree::Vector3 pointPosition, pointVelocityControl, pointForce, pointForceControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;
};



DynamicalComplementarityConstraint::DynamicalComplementarityConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, const std::string &footName, size_t contactIndex, double dissipationGain)
    : iDynTree::optimalcontrol::Constraint (1, "DynamicalComplementarity" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->dissipationGain = dissipationGain;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->velocityControlRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityControlRange.isValid());

    m_pimpl->forceControlRange = controlVariables.getIndexRange(footName + "ForceControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(1);

    m_isLowerBounded = true;
    m_isUpperBounded = true;
    m_lowerBound.zero();
    m_upperBound.zero();

    m_pimpl->stateSparsity.clear();
    m_pimpl->controlSparsity.clear();

    size_t fzIndex = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    size_t pzIndex = static_cast<size_t>(m_pimpl->positionPointRange.offset + 2);

    m_pimpl->stateSparsity.addNonZeroIfNotPresent(0, fzIndex);
    m_pimpl->stateSparsity.addNonZeroIfNotPresent(0, pzIndex);

    size_t fzDotIndex = static_cast<size_t>(m_pimpl->forceControlRange.offset + 2);
    size_t vzIndex = static_cast<size_t>(m_pimpl->velocityControlRange.offset + 2);

    m_pimpl->controlSparsity.addNonZeroIfNotPresent(0, vzIndex);
    m_pimpl->controlSparsity.addNonZeroIfNotPresent(0, fzDotIndex);

}

DynamicalComplementarityConstraint::~DynamicalComplementarityConstraint()
{ }

bool DynamicalComplementarityConstraint::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);

    m_pimpl->constraintValues(0) = m_pimpl->pointForceControl(2) * m_pimpl->pointPosition(2) + m_pimpl->pointForce(2) * m_pimpl->pointVelocityControl(2)
            + m_pimpl->dissipationGain * m_pimpl->pointPosition(2) * m_pimpl->pointForce(2);

    constraint = m_pimpl->constraintValues;

    return true;
}

bool DynamicalComplementarityConstraint::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);

    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);

    m_pimpl->stateJacobianBuffer(0, pzIndex) = m_pimpl->pointForceControl(2) + m_pimpl->dissipationGain * m_pimpl->pointForce(2);
    m_pimpl->stateJacobianBuffer(0, fzIndex) = m_pimpl->pointVelocityControl(2) + m_pimpl->dissipationGain * m_pimpl->pointPosition(2);

    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool DynamicalComplementarityConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);

    unsigned int fzDotIndex = static_cast<unsigned int>(m_pimpl->forceControlRange.offset + 2);
    unsigned int vzIndex = static_cast<unsigned int>(m_pimpl->velocityControlRange.offset + 2);

    m_pimpl->controlJacobianBuffer(0, vzIndex) = m_pimpl->pointForce(2);
    m_pimpl->controlJacobianBuffer(0, fzDotIndex) = m_pimpl->pointPosition(2);

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t DynamicalComplementarityConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t DynamicalComplementarityConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool DynamicalComplementarityConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool DynamicalComplementarityConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}

bool DynamicalComplementarityConstraint::constraintSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                   const iDynTree::VectorDynSize &/*control*/,
                                                                                   const iDynTree::VectorDynSize &lambda,
                                                                                   iDynTree::MatrixDynSize &hessian)
{
    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);

    hessian(fzIndex, pzIndex) = lambda(0) * m_pimpl->dissipationGain;
    hessian(pzIndex, fzIndex) = lambda(0) * m_pimpl->dissipationGain;
    return true;
}

bool DynamicalComplementarityConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                     const iDynTree::VectorDynSize &/*control*/,
                                                                                     const iDynTree::VectorDynSize &/*lambda*/,
                                                                                     iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool DynamicalComplementarityConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/,
                                                                                          const iDynTree::VectorDynSize &/*state*/,
                                                                                          const iDynTree::VectorDynSize &/*control*/,
                                                                                          const iDynTree::VectorDynSize &lambda,
                                                                                          iDynTree::MatrixDynSize &hessian)
{
    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);
    unsigned int fzDotIndex = static_cast<unsigned int>(m_pimpl->forceControlRange.offset + 2);
    unsigned int vzIndex = static_cast<unsigned int>(m_pimpl->velocityControlRange.offset + 2);

    hessian(fzIndex, vzIndex) = lambda(0);
    hessian(pzIndex, fzDotIndex) = lambda(0);
    return true;

}
