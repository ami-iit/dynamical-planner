/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/ContactVelocityControlConstraints.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/MatrixFixSize.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class ContactVelocityControlConstraints::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    HyperbolicSecant activationXY, activationZ;
    iDynTree::Vector3 maximumDerivatives;
    iDynTree::Matrix3x3 dissipationRatios;

    iDynTree::IndexRange positionPointRange, velocityPointRange, forcePointRange, velocityControlRange;
    iDynTree::Vector3 pointPosition, pointVelocity, pointForce, pointVelocityControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;
};


ContactVelocityControlConstraints::ContactVelocityControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                                     const std::string &footName, size_t contactIndex, const HyperbolicSecant &activationXY,
                                                                     const HyperbolicSecant &activationZ,
                                                                     const iDynTree::Vector3 &maximumDerivatives, const iDynTree::Vector3 &dissipationRatios)
    : iDynTree::optimalcontrol::Constraint (6, "VelocityControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->activationXY = activationXY;
    m_pimpl->activationZ = activationZ;
    m_pimpl->maximumDerivatives = maximumDerivatives;
    iDynTree::toEigen(m_pimpl->dissipationRatios) = iDynTree::toEigen(dissipationRatios).asDiagonal();

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->velocityPointRange = stateVariables.getIndexRange(footName + "VelocityPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->velocityControlRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityControlRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(6, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(6, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(6);

    m_isLowerBounded = true;
    m_isUpperBounded = false;
}

bool ContactVelocityControlConstraints::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointVelocity = m_pimpl->stateVariables(m_pimpl->velocityPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaXY = m_pimpl->activationXY.eval(m_pimpl->pointPosition(2));
    double deltaZ = m_pimpl->activationZ.eval(m_pimpl->pointForce(2));

    iDynTree::iDynTreeEigenVector constraintMap = iDynTree::toEigen(m_pimpl->constraintValues);

    constraintMap.topRows<2>() = deltaXY * iDynTree::toEigen(m_pimpl->maximumDerivatives).topRows<2>() - (1-deltaXY) * iDynTree::toEigen(m_pimpl->dissipationRatios).topLeftCorner<2,2>() *
            iDynTree::toEigen(m_pimpl->pointVelocity).topRows<2>() - iDynTree::toEigen(m_pimpl->pointVelocityControl).topRows<2>();

    constraintMap(2) = deltaZ * m_pimpl->maximumDerivatives(2) - (1 - deltaZ) * m_pimpl->dissipationRatios(2,2) * m_pimpl->pointVelocity(2) - m_pimpl->pointVelocityControl(2);

    constraintMap.segment<2>(3) = iDynTree::toEigen(m_pimpl->pointVelocityControl).topRows<2>() + deltaXY * iDynTree::toEigen(m_pimpl->maximumDerivatives).topRows<2>() +
            (1 - deltaXY) * iDynTree::toEigen(m_pimpl->dissipationRatios).topLeftCorner<2,2>() * iDynTree::toEigen(m_pimpl->pointVelocity).topRows<2>();

    constraintMap(5) = m_pimpl->pointVelocityControl(2) + deltaZ * m_pimpl->maximumDerivatives(2) + (1 - deltaZ) * m_pimpl->dissipationRatios(2,2) * m_pimpl->pointVelocity(2);

    constraint = m_pimpl->constraintValues;

    return true;
}

bool ContactVelocityControlConstraints::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointVelocity = m_pimpl->stateVariables(m_pimpl->velocityPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointVelocityControl = m_pimpl->controlVariables(m_pimpl->velocityControlRange);
    double deltaXY = m_pimpl->activationXY.eval(m_pimpl->pointPosition(2));
    double deltaXYDerivative = m_pimpl->activationXY.evalDerivative(m_pimpl->pointPosition(2));
    double deltaZ = m_pimpl->activationZ.eval(m_pimpl->pointForce(2));
    double deltaZDerivative =  m_pimpl->activationZ.evalDerivative(m_pimpl->pointForce(2));


    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);

    jacobianMap.block<2, 1>(0, m_pimpl->positionPointRange.offset + 2) = deltaXYDerivative * iDynTree::toEigen(m_pimpl->maximumDerivatives).topRows<2>() +
            deltaXYDerivative * iDynTree::toEigen(m_pimpl->dissipationRatios).topLeftCorner<2,2>() * iDynTree::toEigen(m_pimpl->pointVelocity).topRows<2>();

    jacobianMap.block<2, 2>(0, m_pimpl->velocityPointRange.offset) = - (1 - deltaXY) * iDynTree::toEigen(m_pimpl->dissipationRatios).topLeftCorner<2,2>();

    jacobianMap(2, m_pimpl->forcePointRange.offset + 2) = deltaZDerivative * m_pimpl->maximumDerivatives(2) + deltaZDerivative * m_pimpl->dissipationRatios(2,2) * m_pimpl->pointVelocity(2);

    jacobianMap(2, m_pimpl->velocityPointRange.offset + 2) = - (1 - deltaZ) * m_pimpl->dissipationRatios(2,2);

    jacobianMap.block<2, 1>(3, m_pimpl->positionPointRange.offset + 2) = deltaXYDerivative * iDynTree::toEigen(m_pimpl->maximumDerivatives).topRows<2>() -
            deltaXYDerivative * iDynTree::toEigen(m_pimpl->dissipationRatios).topLeftCorner<2,2>() * iDynTree::toEigen(m_pimpl->pointVelocity).topRows<2>();

    jacobianMap.block<2, 2>(3, m_pimpl->velocityPointRange.offset) = (1 - deltaXY) * iDynTree::toEigen(m_pimpl->dissipationRatios).topLeftCorner<2,2>();

    jacobianMap(5, m_pimpl->forcePointRange.offset + 2) = deltaZDerivative * m_pimpl->maximumDerivatives(2) - deltaZDerivative * m_pimpl->dissipationRatios(2,2) * m_pimpl->pointVelocity(2);

    jacobianMap(5, m_pimpl->velocityPointRange.offset + 2) = (1 - deltaZ) * m_pimpl->dissipationRatios(2,2);


    jacobian = m_pimpl->stateJacobianBuffer;

    return true;

}

bool ContactVelocityControlConstraints::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<3,3>(0, m_pimpl->velocityControlRange.offset).setIdentity();
    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<3,3>(0, m_pimpl->velocityControlRange.offset) *= -1;

    iDynTree::toEigen(m_pimpl->controlJacobianBuffer).block<3,3>(3, m_pimpl->velocityControlRange.offset).setIdentity();

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ContactVelocityControlConstraints::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ContactVelocityControlConstraints::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}
