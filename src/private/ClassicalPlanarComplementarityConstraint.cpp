/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <DynamicalPlannerPrivate/Constraints/ClassicalPlanarComplementarityConstraint.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>

using namespace DynamicalPlanner::Private;

class ClassicalPlanarComplementarityConstraint::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    double tolerance;

    iDynTree::IndexRange forcePointRange, velocityPointRange;
    iDynTree::Vector3 pointForce, pointVelocity;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateJacobianSparsity, controlJacobianSparsity;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};



ClassicalPlanarComplementarityConstraint::ClassicalPlanarComplementarityConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, const std::string &footName, size_t contactIndex, double tolerance)
    : iDynTree::optimalcontrol::Constraint (2, "ClassicalPlanarComplementarity" + footName + std::to_string(contactIndex))
      , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->tolerance = tolerance;

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->velocityPointRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityPointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(2, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(2, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(2);

    m_isLowerBounded = true;
    m_isUpperBounded = true;
    m_upperBound(0) = m_pimpl->tolerance;
    m_upperBound(1) = m_pimpl->tolerance;
    m_lowerBound(0) = -m_pimpl->tolerance;
    m_lowerBound(1) = -m_pimpl->tolerance;

    m_pimpl->stateJacobianSparsity.clear();
    m_pimpl->controlJacobianSparsity.clear();

    size_t fzIndex = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    size_t vxIndex = static_cast<size_t>(m_pimpl->velocityPointRange.offset);
    size_t vyIndex = static_cast<size_t>(m_pimpl->velocityPointRange.offset + 1);

    m_pimpl->stateJacobianSparsity.add(0, fzIndex);
    m_pimpl->stateJacobianSparsity.add(1, fzIndex);

    m_pimpl->controlJacobianSparsity.add(0, vxIndex);
    m_pimpl->controlJacobianSparsity.add(1, vyIndex);

    m_pimpl->stateHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();
    m_pimpl->mixedHessianSparsity.clear();
    m_pimpl->mixedHessianSparsity.add(fzIndex, vxIndex);
    m_pimpl->mixedHessianSparsity.add(fzIndex, vyIndex);
}

ClassicalPlanarComplementarityConstraint::~ClassicalPlanarComplementarityConstraint()
{ }

bool ClassicalPlanarComplementarityConstraint::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocity = m_pimpl->controlVariables(m_pimpl->velocityPointRange);

    m_pimpl->constraintValues(0) = m_pimpl->pointVelocity(0) * m_pimpl->pointForce(2);
    m_pimpl->constraintValues(1) = m_pimpl->pointVelocity(1) * m_pimpl->pointForce(2);

    constraint = m_pimpl->constraintValues;

    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocity = m_pimpl->controlVariables(m_pimpl->velocityPointRange);

    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    m_pimpl->stateJacobianBuffer(0, fzIndex) = m_pimpl->pointVelocity(0);
    m_pimpl->stateJacobianBuffer(1, fzIndex) = m_pimpl->pointVelocity(1);

    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocity = m_pimpl->controlVariables(m_pimpl->velocityPointRange);

    size_t vxIndex = static_cast<size_t>(m_pimpl->velocityPointRange.offset);
    size_t vyIndex = static_cast<size_t>(m_pimpl->velocityPointRange.offset + 1);

    m_pimpl->controlJacobianBuffer(0, vxIndex) = m_pimpl->pointForce(2);
    m_pimpl->controlJacobianBuffer(1, vyIndex) = m_pimpl->pointForce(2);

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ClassicalPlanarComplementarityConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ClassicalPlanarComplementarityConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool ClassicalPlanarComplementarityConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateJacobianSparsity;
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlJacobianSparsity;
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                   const iDynTree::VectorDynSize &/*control*/,
                                                                                   const iDynTree::VectorDynSize &/*lambda*/,
                                                                                   iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                     const iDynTree::VectorDynSize &/*control*/,
                                                                                     const iDynTree::VectorDynSize &/*lambda*/,
                                                                                     iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/,
                                                                                          const iDynTree::VectorDynSize &/*state*/,
                                                                                          const iDynTree::VectorDynSize &/*control*/,
                                                                                          const iDynTree::VectorDynSize &lambda,
                                                                                          iDynTree::MatrixDynSize &hessian)
{

    size_t fzIndex = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    size_t vxIndex = static_cast<size_t>(m_pimpl->velocityPointRange.offset);
    size_t vyIndex = static_cast<size_t>(m_pimpl->velocityPointRange.offset + 1);

    hessian(fzIndex, vxIndex) = lambda(0);
    hessian(fzIndex, vyIndex) = lambda(1);

    return true;

}

bool ClassicalPlanarComplementarityConstraint::constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ClassicalPlanarComplementarityConstraint::constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
