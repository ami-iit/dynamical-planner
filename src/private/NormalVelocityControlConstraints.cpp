/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Constraints/NormalVelocityControlConstraints.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/MatrixFixSize.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class NormalVelocityControlConstraints::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    HyperbolicTangent normalVelocityActivation;
    double maximumDerivative;

    iDynTree::IndexRange positionPointRange, velocityControlRange;
    iDynTree::Vector3 pointPosition, pointVelocityControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateJacobianSparsity, controlJacobianSparsity;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

};


NormalVelocityControlConstraints::NormalVelocityControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                                   const std::string &footName, size_t contactIndex, const HyperbolicTangent &lowerBoundMultiplier, double maximumDerivative)
    : iDynTree::optimalcontrol::Constraint (2, "NormalVelocityControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->normalVelocityActivation = lowerBoundMultiplier;
    m_pimpl->maximumDerivative = maximumDerivative;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->velocityControlRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityControlRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(2, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(2, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(2);

    m_isLowerBounded = true;
    m_isUpperBounded = false;
    m_lowerBound.zero();

    m_pimpl->stateJacobianSparsity.clear();
    m_pimpl->controlJacobianSparsity.clear();

    size_t pzIndex = static_cast<size_t>(m_pimpl->positionPointRange.offset + 2);
    m_pimpl->stateJacobianSparsity.addDenseBlock(0, pzIndex, 2, 1);
    size_t velocityIndex = static_cast<size_t>(m_pimpl->velocityControlRange.offset + 2);
    m_pimpl->controlJacobianSparsity.addDenseBlock(0, velocityIndex, 2, 1);

    m_pimpl->stateHessianSparsity.add(pzIndex, pzIndex);
    m_pimpl->mixedHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();
}

NormalVelocityControlConstraints::~NormalVelocityControlConstraints()
{ }

bool NormalVelocityControlConstraints::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaZ = m_pimpl->normalVelocityActivation.eval(m_pimpl->pointPosition(2));

    iDynTree::iDynTreeEigenVector constraintMap = iDynTree::toEigen(m_pimpl->constraintValues);

    constraintMap(0) = m_pimpl->maximumDerivative - m_pimpl->pointVelocityControl(2);

    constraintMap(1) = m_pimpl->pointVelocityControl(2) + deltaZ * m_pimpl->maximumDerivative;

    constraint = m_pimpl->constraintValues;

    return true;
}

bool NormalVelocityControlConstraints::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaZDerivative =  m_pimpl->normalVelocityActivation.evalDerivative(m_pimpl->pointPosition(2));

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);

    m_pimpl->stateJacobianBuffer(1, pzIndex) = deltaZDerivative * m_pimpl->maximumDerivative;

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;

}

bool NormalVelocityControlConstraints::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    unsigned int velocityIndex = static_cast<unsigned int>(m_pimpl->velocityControlRange.offset + 2);

    m_pimpl->controlJacobianBuffer(0, velocityIndex) = -1.0;

    m_pimpl->controlJacobianBuffer(1, velocityIndex) = 1.0;

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t NormalVelocityControlConstraints::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t NormalVelocityControlConstraints::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool NormalVelocityControlConstraints::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateJacobianSparsity;
    return true;
}

bool NormalVelocityControlConstraints::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlJacobianSparsity;
    return true;
}

bool NormalVelocityControlConstraints::constraintSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state,
                                                                                 const iDynTree::VectorDynSize &/*control*/,
                                                                                 const iDynTree::VectorDynSize &lambda,
                                                                                 iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);
    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    double deltaDoubleDerivative = m_pimpl->normalVelocityActivation.evalDoubleDerivative(m_pimpl->pointPosition(2));

    hessian(pzIndex, pzIndex) = (lambda(1)) * deltaDoubleDerivative * m_pimpl->maximumDerivative;

    return true;
}

bool NormalVelocityControlConstraints::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                   const iDynTree::VectorDynSize &/*control*/,
                                                                                   const iDynTree::VectorDynSize &/*lambda*/,
                                                                                   iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool NormalVelocityControlConstraints::constraintSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                        const iDynTree::VectorDynSize &/*control*/,
                                                                                        const iDynTree::VectorDynSize &/*lambda*/,
                                                                                        iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool NormalVelocityControlConstraints::constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool NormalVelocityControlConstraints::constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool NormalVelocityControlConstraints::constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
