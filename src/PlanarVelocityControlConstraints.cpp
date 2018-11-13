/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/PlanarVelocityControlConstraints.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/MatrixFixSize.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class PlanarVelocityControlConstraints::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    HyperbolicTangent planarVelocityActivation;
    iDynTree::Vector2 maximumDerivatives;

    iDynTree::IndexRange positionPointRange, forcePointRange, velocityControlRange;
    iDynTree::Vector3 pointPosition, pointForce, pointVelocityControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;
};


PlanarVelocityControlConstraints::PlanarVelocityControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                                   const std::string &footName, size_t contactIndex,
                                                                   const HyperbolicTangent &planarVelocityActivation, double xMaximumDerivative, double yMaximumDerivative)
    : iDynTree::optimalcontrol::Constraint (4, "PlanarVelocityControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->planarVelocityActivation = planarVelocityActivation;
    m_pimpl->maximumDerivatives(0) = xMaximumDerivative;
    m_pimpl->maximumDerivatives(1) = yMaximumDerivative;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->velocityControlRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityControlRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(4, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(4, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(4);

    m_isLowerBounded = true;
    m_isUpperBounded = false;
    m_lowerBound.zero();

    m_pimpl->stateSparsity.clear();
    m_pimpl->controlSparsity.clear();

    size_t pzIndex = static_cast<size_t>(m_pimpl->positionPointRange.offset + 2);

    m_pimpl->stateSparsity.addDenseBlock(0, pzIndex, 4, 1);

    m_pimpl->controlSparsity.addIdentityBlock(0, static_cast<size_t>(m_pimpl->velocityControlRange.offset), 2);
    m_pimpl->controlSparsity.addIdentityBlock(2, static_cast<size_t>(m_pimpl->velocityControlRange.offset), 2);
}

PlanarVelocityControlConstraints::~PlanarVelocityControlConstraints()
{ }

bool PlanarVelocityControlConstraints::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaXY = m_pimpl->planarVelocityActivation.eval(m_pimpl->pointPosition(2));

    iDynTree::iDynTreeEigenVector constraintMap = iDynTree::toEigen(m_pimpl->constraintValues);

    constraintMap.topRows<2>() = deltaXY * iDynTree::toEigen(m_pimpl->maximumDerivatives) - iDynTree::toEigen(m_pimpl->pointVelocityControl).topRows<2>();

    constraintMap.bottomRows<2>() = iDynTree::toEigen(m_pimpl->pointVelocityControl).topRows<2>() + deltaXY * iDynTree::toEigen(m_pimpl->maximumDerivatives);

    constraint = m_pimpl->constraintValues;

    return true;
}

bool PlanarVelocityControlConstraints::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaXYDerivative = m_pimpl->planarVelocityActivation.evalDerivative(m_pimpl->pointPosition(2));

    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);

    jacobianMap.block<2, 1>(0, m_pimpl->positionPointRange.offset + 2) = deltaXYDerivative * iDynTree::toEigen(m_pimpl->maximumDerivatives);

    jacobianMap.block<2, 1>(2, m_pimpl->positionPointRange.offset + 2) = deltaXYDerivative * iDynTree::toEigen(m_pimpl->maximumDerivatives);

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;

}

bool PlanarVelocityControlConstraints::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<2,2>(0, m_pimpl->velocityControlRange.offset).setIdentity();
    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<2,2>(0, m_pimpl->velocityControlRange.offset) *= -1;

    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<2,2>(2, m_pimpl->velocityControlRange.offset).setIdentity();

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t PlanarVelocityControlConstraints::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t PlanarVelocityControlConstraints::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool PlanarVelocityControlConstraints::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool PlanarVelocityControlConstraints::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}
