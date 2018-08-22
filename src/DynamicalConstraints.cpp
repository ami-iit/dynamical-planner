/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/DynamicalConstraints.h>
#include <private/SpanUtils.h>
#include <private/QuaternionUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/Utils.h>
#include <cassert>
#include <iostream>
#include <regex>

using namespace DynamicalPlanner::Private;

class DynamicalConstraints::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;
    VariablesLabeller dynamics;
    double totalMass;
    size_t leftPoints, rightPoints;
    iDynTree::Vector6 gravityVector;
    iDynTree::Rotation baseRotation;

    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    void computeFootRelatedDynamics(const std::string& footName, size_t numberOfPoints) {
        Eigen::Vector3d distance, appliedForce;
        std::string positionLabel, velocityLabel, forceLabel, velocityControlLabel, forceControlLabel;

        for (size_t i; i < numberOfPoints; ++i) {
            positionLabel = footName + "PositionPoint" + std::to_string(i);
            velocityLabel = footName + "VelocityPoint" + std::to_string(i);
            forceLabel = footName + "ForcePoint" + std::to_string(i);
            velocityControlLabel = footName + "VelocityControlPoint" + std::to_string(i);
            forceControlLabel = footName + "ForceControlPoint" + std::to_string(i);

            spanToEigen(dynamics(forceLabel)) = spanToEigen(controlVariables(forceControlLabel));
            spanToEigen(dynamics(velocityLabel)) = spanToEigen(controlVariables(velocityControlLabel));
            spanToEigen(dynamics(positionLabel)) = spanToEigen(stateVariables(velocityLabel));

            spanToEigen(dynamics("Momentum")).topRows<3>() +=  spanToEigen(stateVariables(forceLabel));

            distance = spanToEigen(stateVariables(positionLabel)) - spanToEigen(stateVariables("CoMPosition"));
            appliedForce = spanToEigen(stateVariables(forceLabel));
            spanToEigen(dynamics("Momentum")).bottomRows<3>() += distance.cross(appliedForce);
        }
    }

    void computeFootRelatedStateJacobian(const std::string& footName, size_t numberOfPoints) {
        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(stateJacobianBuffer);
        iDynTree::IndexRange rowRange, columnRange;
        Eigen::Vector3d distance, appliedForce;
        std::string positionLabel, velocityLabel, forceLabel;

        for (size_t i; i < numberOfPoints; ++i) {
            positionLabel = footName + "PositionPoint" + std::to_string(i);
            velocityLabel = footName + "VelocityPoint" + std::to_string(i);
            forceLabel = footName + "ForcePoint" + std::to_string(i);

            rowRange = stateVariables.getIndexRange(positionLabel);
            columnRange = stateVariables.getIndexRange(velocityLabel);
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).setIdentity();

            rowRange = stateVariables.getIndexRange("Momentum");
            columnRange = stateVariables.getIndexRange(forceLabel);
            distance = spanToEigen(stateVariables(positionLabel)) - spanToEigen(stateVariables("CoMPosition"));
            appliedForce = spanToEigen(stateVariables(forceLabel));
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).topRows<3>().setIdentity();
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).bottomRows<3>() = iDynTree::skew(distance);
            columnRange = stateVariables.getIndexRange(positionLabel);
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).bottomRows<3>() = -iDynTree::skew(appliedForce);
            columnRange = stateVariables.getIndexRange("CoMPosition");
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).bottomRows<3>() = iDynTree::skew(appliedForce);
        }
    }
};


DynamicalConstraints::DynamicalConstraints(VariablesLabeller &stateVariables, VariablesLabeller &controlVariables, double totalMass)
   : iDynTree::optimalcontrol::DynamicalSystem (stateVariables.size(), controlVariables.size())
   , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;
    m_pimpl->dynamics = m_pimpl->stateVariables;
    m_pimpl->totalMass = totalMass;

    m_pimpl->gravityVector.zero();
    m_pimpl->gravityVector(2) = -9.81;

    m_pimpl->leftPoints = 0;
    m_pimpl->rightPoints = 0;

    m_pimpl->stateJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    for (auto label : stateVariables.listOfLabels()) {
        if (label.find("LeftForcePoint") != std::string::npos) {
            m_pimpl->leftPoints++;
        }

        if (label.find("RightForcePoint") != std::string::npos) {
             m_pimpl->rightPoints++;
        }
    }

    //Check that all corresponding variables are available
    for (size_t i; i < m_pimpl->leftPoints; ++i) {
        if (!stateVariables.getIndexRange("LeftForcePoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftForcePoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("LeftVelocityPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftVelocityPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("LeftPositionPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftPositionPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("LeftForceControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftForceControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("LeftVelocityControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftVelocityControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }
    }

    for (size_t i; i < m_pimpl->rightPoints; ++i) {
        if (!stateVariables.getIndexRange("RightForcePoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightForcePoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("RightVelocityPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightVelocityPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("RightPositionPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightPositionPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("RightForceControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightForceControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("RightVelocityControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightVelocityControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }
    }

}

DynamicalConstraints::~DynamicalConstraints()
{

}

bool DynamicalConstraints::dynamics(const iDynTree::VectorDynSize &state, double /*time*/, iDynTree::VectorDynSize &stateDynamics)
{
    spanToEigen(m_pimpl->stateVariables.values()) = iDynTree::toEigen(state); //this line must remain before those computing the feet related quantities
    spanToEigen(m_pimpl->dynamics.values()) = iDynTree::toEigen(state); //this line must remain before those computing the feet related quantities
    spanToEigen(m_pimpl->controlVariables.values()) = iDynTree::toEigen(controlInput()); //this line must remain before those computing the feet related quantities


    spanToEigen(m_pimpl->dynamics("Momentum")) = m_pimpl->totalMass * iDynTree::toEigen(m_pimpl->gravityVector); //this line must remain before those computing the feet related quantities

    m_pimpl->computeFootRelatedDynamics("Left", m_pimpl->leftPoints);

    m_pimpl->computeFootRelatedDynamics("Right", m_pimpl->rightPoints);

    spanToEigen(m_pimpl->dynamics("CoMPosition")) = spanToEigen(m_pimpl->stateVariables("Momentum")).topRows<3>()/m_pimpl->totalMass;

    iDynTree::Vector4 normalizedQuaternion;
    iDynTree::toEigen(normalizedQuaternion) = spanToEigen(m_pimpl->stateVariables("BaseQuaternion")).normalized();

    assert(QuaternionBoundsRespected(normalizedQuaternion));

    m_pimpl->baseRotation.fromQuaternion(normalizedQuaternion);

    spanToEigen(m_pimpl->dynamics("BasePosition")) = iDynTree::toEigen(m_pimpl->baseRotation) * spanToEigen(m_pimpl->controlVariables("BaseVelocity")).topRows<3>();

    spanToEigen(m_pimpl->dynamics("BaseQuaternion")) = iDynTree::toEigen(QuaternionLeftTrivializedDerivative(normalizedQuaternion)) * spanToEigen(m_pimpl->controlVariables("BaseVelocity")).bottomRows<3>();

    spanToEigen(m_pimpl->dynamics("JointsPosition")) = spanToEigen(m_pimpl->controlVariables("JointsVelocity"));

    iDynTree::toEigen(stateDynamics) = spanToEigen(m_pimpl->dynamics.values());

    return true;
}

bool DynamicalConstraints::dynamicsStateFirstDerivative(const iDynTree::VectorDynSize &state, double /*time*/, iDynTree::MatrixDynSize &dynamicsDerivative)
{
    spanToEigen(m_pimpl->stateVariables.values()) = iDynTree::toEigen(state);
    spanToEigen(m_pimpl->controlVariables.values()) = iDynTree::toEigen(controlInput());
    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);

    m_pimpl->computeFootRelatedStateJacobian("Left", m_pimpl->leftPoints);
    m_pimpl->computeFootRelatedStateJacobian("Right", m_pimpl->rightPoints);

    iDynTree::IndexRange rowRange, columnRange;
    rowRange = m_pimpl->stateVariables.getIndexRange("CoMPosition");
    columnRange = m_pimpl->stateVariables.getIndexRange("Momentum");
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).leftCols<3>().setIdentity();
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).leftCols<3>() *= 1/m_pimpl->totalMass;

    rowRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    columnRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");

    iDynTree::Vector4 normalizedQuaternion, quaternion;
    iDynTree::toEigen(quaternion) = spanToEigen(m_pimpl->stateVariables("BaseQuaternion"));
    normalizedQuaternion = NormailizedQuaternion(quaternion);
    iDynTree::Vector3 baseLinearVelocity;
    iDynTree::toEigen(baseLinearVelocity) = spanToEigen(m_pimpl->controlVariables("BaseVelocity")).topRows<3>();
    iDynTree::Matrix4x4 normalizedQuaternionDerivative = NormalizedQuaternionDerivative(quaternion);
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size) = iDynTree::toEigen(RotatedVectorQuaternionJacobian(baseLinearVelocity, normalizedQuaternion)) * iDynTree::toEigen(normalizedQuaternionDerivative);

    rowRange = columnRange;
    iDynTree::Vector3 baseAngularVelocity;
    iDynTree::toEigen(baseAngularVelocity) = spanToEigen(m_pimpl->controlVariables("BaseVelocity")).bottomRows<3>();
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeTimesOmegaJacobian(baseAngularVelocity)) * iDynTree::toEigen(normalizedQuaternionDerivative);

    dynamicsDerivative = m_pimpl->stateJacobianBuffer;
    return true;
}

bool DynamicalConstraints::dynamicsControlFirstDerivative(const iDynTree::VectorDynSize &state, double time, iDynTree::MatrixDynSize &dynamicsDerivative)
{

}
