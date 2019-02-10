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
    HyperbolicSecant normalVelocityActivation;
    double maximumDerivative;

    iDynTree::IndexRange positionPointRange, forcePointRange, velocityControlRange;
    iDynTree::Vector3 pointPosition, pointForce, pointVelocityControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;
};


NormalVelocityControlConstraints::NormalVelocityControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                                   const std::string &footName, size_t contactIndex, const HyperbolicSecant &normalVelocityActivation, double maximumDerivative)
    : iDynTree::optimalcontrol::Constraint (2, "NormalVelocityControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->normalVelocityActivation = normalVelocityActivation;
    m_pimpl->maximumDerivative = maximumDerivative;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

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

    m_pimpl->stateSparsity.clear();
    m_pimpl->controlSparsity.clear();

    size_t forceIndex = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    m_pimpl->stateSparsity.addDenseBlock(0, forceIndex, 2, 1);
    size_t velocityIndex = static_cast<size_t>(m_pimpl->velocityControlRange.offset + 2);
    m_pimpl->controlSparsity.addDenseBlock(0, velocityIndex, 2, 1);
}

NormalVelocityControlConstraints::~NormalVelocityControlConstraints()
{ }

bool NormalVelocityControlConstraints::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaZ = m_pimpl->normalVelocityActivation.eval(m_pimpl->pointForce(2));

    iDynTree::iDynTreeEigenVector constraintMap = iDynTree::toEigen(m_pimpl->constraintValues);

    constraintMap(0) = deltaZ * m_pimpl->maximumDerivative - m_pimpl->pointVelocityControl(2);

    constraintMap(1) = m_pimpl->pointVelocityControl(2) + deltaZ * m_pimpl->maximumDerivative;

    constraint = m_pimpl->constraintValues;

    return true;
}

bool NormalVelocityControlConstraints::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaZDerivative =  m_pimpl->normalVelocityActivation.evalDerivative(m_pimpl->pointForce(2));

    unsigned int forceIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    m_pimpl->stateJacobianBuffer(0, forceIndex) = deltaZDerivative * m_pimpl->maximumDerivative;

    m_pimpl->stateJacobianBuffer(1, forceIndex) = deltaZDerivative * m_pimpl->maximumDerivative;

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
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool NormalVelocityControlConstraints::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}
