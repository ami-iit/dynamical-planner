/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <DynamicalPlannerPrivate/Constraints/ClassicalComplementarityConstraint.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>

using namespace DynamicalPlanner::Private;

class ClassicalComplementarityConstraint::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    double tolerance;

    iDynTree::IndexRange positionPointRange, forcePointRange;
    iDynTree::Vector3 pointPosition, pointForce;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateJacobianSparsity, controlJacobianSparsity;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};



ClassicalComplementarityConstraint::ClassicalComplementarityConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, const std::string &footName, size_t contactIndex, double tolerance)
    : iDynTree::optimalcontrol::Constraint (1, "ClassicalComplementarity" + footName + std::to_string(contactIndex))
      , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->tolerance = tolerance;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(1);

    m_isLowerBounded = false;
    m_isUpperBounded = true;
    m_lowerBound.zero();
    m_upperBound(0) = m_pimpl->tolerance;

    m_pimpl->stateJacobianSparsity.clear();
    m_pimpl->controlJacobianSparsity.clear();

    size_t fzIndex = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    size_t pzIndex = static_cast<size_t>(m_pimpl->positionPointRange.offset + 2);

    m_pimpl->stateJacobianSparsity.add(0, fzIndex);
    m_pimpl->stateJacobianSparsity.add(0, pzIndex);

    m_pimpl->stateHessianSparsity.add(fzIndex, pzIndex);
    m_pimpl->stateHessianSparsity.add(pzIndex, fzIndex);

    m_pimpl->controlHessianSparsity.clear();
    m_pimpl->mixedHessianSparsity.clear();
}

ClassicalComplementarityConstraint::~ClassicalComplementarityConstraint()
{ }

bool ClassicalComplementarityConstraint::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);


    m_pimpl->constraintValues(0) = m_pimpl->pointPosition(2) * m_pimpl->pointForce(2);

    constraint = m_pimpl->constraintValues;

    return true;
}

bool ClassicalComplementarityConstraint::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);

    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);

    m_pimpl->stateJacobianBuffer(0, pzIndex) = m_pimpl->pointForce(2);
    m_pimpl->stateJacobianBuffer(0, fzIndex) = m_pimpl->pointPosition(2);

    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool ClassicalComplementarityConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ClassicalComplementarityConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ClassicalComplementarityConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool ClassicalComplementarityConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateJacobianSparsity;
    return true;
}

bool ClassicalComplementarityConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlJacobianSparsity;
    return true;
}

bool ClassicalComplementarityConstraint::constraintSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                   const iDynTree::VectorDynSize &/*control*/,
                                                                                   const iDynTree::VectorDynSize &lambda,
                                                                                   iDynTree::MatrixDynSize &hessian)
{
    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);

    hessian(fzIndex, pzIndex) = lambda(0);
    hessian(pzIndex, fzIndex) = lambda(0);
    return true;
}

bool ClassicalComplementarityConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                     const iDynTree::VectorDynSize &/*control*/,
                                                                                     const iDynTree::VectorDynSize &/*lambda*/,
                                                                                     iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ClassicalComplementarityConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/,
                                                                                          const iDynTree::VectorDynSize &/*state*/,
                                                                                          const iDynTree::VectorDynSize &/*control*/,
                                                                                          const iDynTree::VectorDynSize &/*lambda*/,
                                                                                          iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;

}

bool ClassicalComplementarityConstraint::constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ClassicalComplementarityConstraint::constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ClassicalComplementarityConstraint::constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
