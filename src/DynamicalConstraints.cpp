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
    //size_t leftPoints, rightPoints;
    iDynTree::Vector6 gravityVector;
    iDynTree::Rotation baseRotation;

    typedef struct {
        std::vector<iDynTree::IndexRange> positionPoints, velocityPoints, forcePoints, velocityControlPoints, forceControlPoints;
    } FootRanges;
    FootRanges leftRanges, rightRanges;
    iDynTree::IndexRange momentumRange, comPositionRange, basePositionRange, baseQuaternionRange, jointsPositionRange, baseVelocityRange, jointsVelocityRange;

    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    void checkFootVariables(const std::string& footName, size_t numberOfPoints, FootRanges& foot) {
        foot.positionPoints.resize(numberOfPoints);
        foot.velocityPoints.resize(numberOfPoints);
        foot.forcePoints.resize(numberOfPoints);
        foot.velocityControlPoints.resize(numberOfPoints);
        foot.forceControlPoints.resize(numberOfPoints);

        iDynTree::IndexRange obtainedRange;
        for (size_t i; i < numberOfPoints; ++i) {
            obtainedRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(i));

            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "ForcePoint" + std::to_string(i) << " not available among state variables." << std::endl;
                assert(false);
            }
            foot.forcePoints[i] = obtainedRange;

            obtainedRange = stateVariables.getIndexRange(footName + "VelocityPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "VelocityPoint" + std::to_string(i) << " not available among state variables." << std::endl;
                assert(false);
            }
            foot.velocityPoints[i] = obtainedRange;

            obtainedRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "PositionPoint" + std::to_string(i) << " not available among state variables." << std::endl;
                assert(false);
            }
            foot.positionPoints[i] = obtainedRange;

            obtainedRange = controlVariables.getIndexRange(footName + "ForceControlPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "ForceControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
                assert(false);
            }
            foot.forceControlPoints[i] = obtainedRange;

            obtainedRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "VelocityControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
                assert(false);
            }
            foot.velocityControlPoints[i] = obtainedRange;
        }
    }


    void computeFootRelatedDynamics(FootRanges& foot) {
        Eigen::Vector3d distance, appliedForce;

        for (size_t i; i < foot.positionPoints.size(); ++i) {

            spanToEigen(dynamics(foot.forcePoints[i])) = spanToEigen(controlVariables(foot.forceControlPoints[i]));
            spanToEigen(dynamics(foot.velocityPoints[i])) = spanToEigen(controlVariables(foot.velocityControlPoints[i]));
            spanToEigen(dynamics(foot.positionPoints[i])) = spanToEigen(stateVariables(foot.velocityPoints[i]));

            spanToEigen(dynamics(momentumRange)).topRows<3>() +=  spanToEigen(stateVariables(foot.forcePoints[i]));

            distance = spanToEigen(stateVariables(foot.positionPoints[i])) - spanToEigen(stateVariables(comPositionRange));
            appliedForce = spanToEigen(stateVariables(foot.forcePoints[i]));
            spanToEigen(dynamics(momentumRange)).bottomRows<3>() += distance.cross(appliedForce);
        }
    }

    void computeFootRelatedStateJacobian(FootRanges& foot) {
        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(stateJacobianBuffer);
        iDynTree::IndexRange rowRange, columnRange;
        Eigen::Vector3d distance, appliedForce;


        for (size_t i; i < foot.positionPoints.size(); ++i) {

            rowRange = foot.positionPoints[i];
            columnRange = foot.velocityPoints[i];
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).setIdentity();

            rowRange = momentumRange;
            columnRange = foot.forcePoints[i];
            distance = spanToEigen(stateVariables(foot.positionPoints[i])) - spanToEigen(stateVariables(comPositionRange));
            appliedForce = spanToEigen(stateVariables(foot.forcePoints[i]));
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).topRows<3>().setIdentity();
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).bottomRows<3>() = iDynTree::skew(distance);
            columnRange = foot.positionPoints[i];
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).bottomRows<3>() = -iDynTree::skew(appliedForce);
            columnRange = comPositionRange;
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).bottomRows<3>() = iDynTree::skew(appliedForce);
        }
    }

    void computeFootRelatedControlJacobian(FootRanges& foot) {
        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(controlJacobianBuffer);
        iDynTree::IndexRange rowRange, columnRange;

        for (size_t i; i < foot.positionPoints.size(); ++i) {
            rowRange = foot.forcePoints[i];
            columnRange = foot.forceControlPoints[i];
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).setIdentity();

            rowRange = foot.velocityPoints[i];
            columnRange = foot.velocityControlPoints[i];
            jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).setIdentity();
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

    m_pimpl->stateJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    size_t leftPoints = 0, rightPoints = 0;
    for (auto label : stateVariables.listOfLabels()) {
        if (label.find("LeftForcePoint") != std::string::npos) {
            leftPoints++;
        }

        if (label.find("RightForcePoint") != std::string::npos) {
             rightPoints++;
        }
    }

    m_pimpl->checkFootVariables("Left",leftPoints, m_pimpl->leftRanges);
    m_pimpl->checkFootVariables("Right", rightPoints, m_pimpl->rightRanges);

    m_pimpl->momentumRange = m_pimpl->stateVariables.getIndexRange("Momentum");
    assert(m_pimpl->momentumRange.isValid());

    m_pimpl->comPositionRange = m_pimpl->stateVariables.getIndexRange("CoMPosition");
    assert(m_pimpl->comPositionRange.isValid());

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->baseVelocityRange = m_pimpl->controlVariables.getIndexRange("BaseVelocity");
    assert(m_pimpl->baseVelocityRange.isValid());

    m_pimpl->jointsVelocityRange = m_pimpl->controlVariables.getIndexRange("JointsVelocity");
    assert(m_pimpl->jointsVelocityRange.isValid());

}

DynamicalConstraints::~DynamicalConstraints()
{

}

bool DynamicalConstraints::dynamics(const iDynTree::VectorDynSize &state, double /*time*/, iDynTree::VectorDynSize &stateDynamics)
{
    spanToEigen(m_pimpl->stateVariables.values()) = iDynTree::toEigen(state); //this line must remain before those computing the feet related quantities
    spanToEigen(m_pimpl->dynamics.values()) = iDynTree::toEigen(state); //this line must remain before those computing the feet related quantities
    spanToEigen(m_pimpl->controlVariables.values()) = iDynTree::toEigen(controlInput()); //this line must remain before those computing the feet related quantities


    spanToEigen(m_pimpl->dynamics(m_pimpl->momentumRange)) = m_pimpl->totalMass * iDynTree::toEigen(m_pimpl->gravityVector); //this line must remain before those computing the feet related quantities

    m_pimpl->computeFootRelatedDynamics(m_pimpl->leftRanges);

    m_pimpl->computeFootRelatedDynamics(m_pimpl->rightRanges);

    spanToEigen(m_pimpl->dynamics(m_pimpl->comPositionRange)) = spanToEigen(m_pimpl->stateVariables(m_pimpl->momentumRange)).topRows<3>()/m_pimpl->totalMass;

    iDynTree::Vector4 normalizedQuaternion;
    iDynTree::toEigen(normalizedQuaternion) = spanToEigen(m_pimpl->stateVariables(m_pimpl->baseQuaternionRange)).normalized();

    assert(QuaternionBoundsRespected(normalizedQuaternion));

    m_pimpl->baseRotation.fromQuaternion(normalizedQuaternion);

    spanToEigen(m_pimpl->dynamics(m_pimpl->basePositionRange)) = iDynTree::toEigen(m_pimpl->baseRotation) * spanToEigen(m_pimpl->controlVariables(m_pimpl->baseVelocityRange)).topRows<3>();

    spanToEigen(m_pimpl->dynamics(m_pimpl->baseQuaternionRange)) = iDynTree::toEigen(QuaternionLeftTrivializedDerivative(normalizedQuaternion)) * spanToEigen(m_pimpl->controlVariables(m_pimpl->baseVelocityRange)).bottomRows<3>();

    spanToEigen(m_pimpl->dynamics(m_pimpl->jointsPositionRange)) = spanToEigen(m_pimpl->controlVariables(m_pimpl->jointsVelocityRange));

    iDynTree::toEigen(stateDynamics) = spanToEigen(m_pimpl->dynamics.values());

    return true;
}

bool DynamicalConstraints::dynamicsStateFirstDerivative(const iDynTree::VectorDynSize &state, double /*time*/, iDynTree::MatrixDynSize &dynamicsDerivative)
{
    spanToEigen(m_pimpl->stateVariables.values()) = iDynTree::toEigen(state);
    spanToEigen(m_pimpl->controlVariables.values()) = iDynTree::toEigen(controlInput());
    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);

    m_pimpl->computeFootRelatedStateJacobian(m_pimpl->leftRanges);
    m_pimpl->computeFootRelatedStateJacobian(m_pimpl->rightRanges);

    iDynTree::IndexRange rowRange, columnRange;
    rowRange = m_pimpl->comPositionRange;
    columnRange = m_pimpl->momentumRange;
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).leftCols<3>().setIdentity();
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).leftCols<3>() *= 1/m_pimpl->totalMass;

    rowRange = m_pimpl->basePositionRange;
    columnRange = m_pimpl->baseQuaternionRange;

    iDynTree::Vector4 normalizedQuaternion, quaternion;
    iDynTree::toEigen(quaternion) = spanToEigen(m_pimpl->stateVariables(m_pimpl->baseQuaternionRange));
    normalizedQuaternion = NormailizedQuaternion(quaternion);
    assert(QuaternionBoundsRespected(normalizedQuaternion));

    iDynTree::Vector3 baseLinearVelocity;
    iDynTree::toEigen(baseLinearVelocity) = spanToEigen(m_pimpl->controlVariables(m_pimpl->baseVelocityRange)).topRows<3>();
    iDynTree::Matrix4x4 normalizedQuaternionDerivative = NormalizedQuaternionDerivative(quaternion);
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size) = iDynTree::toEigen(RotatedVectorQuaternionJacobian(baseLinearVelocity, normalizedQuaternion)) * iDynTree::toEigen(normalizedQuaternionDerivative);

    rowRange = columnRange;
    iDynTree::Vector3 baseAngularVelocity;
    iDynTree::toEigen(baseAngularVelocity) = spanToEigen(m_pimpl->controlVariables(m_pimpl->baseVelocityRange)).bottomRows<3>();
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeTimesOmegaJacobian(baseAngularVelocity)) * iDynTree::toEigen(normalizedQuaternionDerivative);

    dynamicsDerivative = m_pimpl->stateJacobianBuffer;
    return true;
}

bool DynamicalConstraints::dynamicsControlFirstDerivative(const iDynTree::VectorDynSize &state, double /*time*/, iDynTree::MatrixDynSize &dynamicsDerivative)
{
    spanToEigen(m_pimpl->stateVariables.values()) = iDynTree::toEigen(state);
    spanToEigen(m_pimpl->controlVariables.values()) = iDynTree::toEigen(controlInput());
    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->controlJacobianBuffer);

    m_pimpl->computeFootRelatedControlJacobian(m_pimpl->leftRanges);
    m_pimpl->computeFootRelatedControlJacobian(m_pimpl->rightRanges);

    iDynTree::IndexRange rowRange, columnRange;

    iDynTree::Vector4 normalizedQuaternion;
    iDynTree::toEigen(normalizedQuaternion) = spanToEigen(m_pimpl->stateVariables(m_pimpl->baseQuaternionRange)).normalized();
    assert(QuaternionBoundsRespected(normalizedQuaternion));

    m_pimpl->baseRotation.fromQuaternion(normalizedQuaternion);
    rowRange = m_pimpl->basePositionRange;
    columnRange = m_pimpl->baseVelocityRange;
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).leftCols<3>() = iDynTree::toEigen(m_pimpl->baseRotation);

    rowRange = m_pimpl->baseQuaternionRange;
    columnRange = m_pimpl->baseVelocityRange;
    jacobianMap.block(rowRange.offset, columnRange.offset, rowRange.size, columnRange.size).rightCols<3>() = iDynTree::toEigen(QuaternionLeftTrivializedDerivative(normalizedQuaternion));


    dynamicsDerivative = m_pimpl->controlJacobianBuffer;
    return true;
}
