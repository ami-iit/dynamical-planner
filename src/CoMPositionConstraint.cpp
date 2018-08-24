/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/CoMPositionConstraint.h>
#include <private/QuaternionUtils.h>

#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/Core/MatrixDynSize.h>

#include <cassert>

using namespace DynamicalPlanner::Private;

class CoMPositionConstraint::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;
    iDynTree::IndexRange comPositionRange, jointsPositionRange, basePositionRange, baseQuaternionRange;
    iDynTree::VectorDynSize constraintValueBuffer;
    iDynTree::Position basePosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::Rotation baseRotation;
    iDynTree::Vector3 gravity;
    iDynTree::MatrixDynSize comJacobianBuffer, stateJacobianBuffer;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputation> sharedKinDyn;

    void updateRobotState() {

        robotState = sharedKinDyn->currentState();

        iDynTree::toEigen(basePosition) = iDynTree::toEigen(stateVariables(basePositionRange));
        baseQuaternion = stateVariables(baseQuaternionRange);
        baseQuaternionNormalized = NormailizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);

        robotState.world_T_base.setRotation(baseRotation);
        robotState.world_T_base.setPosition(basePosition);

        robotState.s = stateVariables(jointsPositionRange);
    }
};


CoMPositionConstraint::CoMPositionConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, std::shared_ptr<SharedKinDynComputation> sharedKinDyn)
    : iDynTree::optimalcontrol::Constraint (3, "CoMPosition")
    , m_pimpl(new Implementation)
{
    assert(sharedKinDyn);
    assert(sharedKinDyn->isValid());
    m_pimpl->sharedKinDyn = sharedKinDyn;

    m_pimpl->gravity.zero();
    m_pimpl->gravity(2) = -9.81;


    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_isLowerBounded = true;
    m_isUpperBounded = true;

    m_pimpl->comPositionRange = m_pimpl->stateVariables.getIndexRange("CoMPosition");
    assert(m_pimpl->comPositionRange.isValid());

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->constraintValueBuffer.resize(3);
    m_pimpl->constraintValueBuffer.zero();
    m_pimpl->comJacobianBuffer.resize(3, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateJacobianBuffer.resize(3, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->robotState = sharedKinDyn->currentState();
}

bool CoMPositionConstraint::evaluateConstraint(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;

    m_pimpl->updateRobotState();

    iDynTree::toEigen(m_pimpl->constraintValueBuffer) = iDynTree::toEigen(m_pimpl->sharedKinDyn->getCenterOfMassPosition(m_pimpl->robotState)) - iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->comPositionRange));

    constraint = m_pimpl->constraintValueBuffer;

    return true;

}

bool CoMPositionConstraint::constraintJacobianWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;

    m_pimpl->updateRobotState();

    bool ok = m_pimpl->sharedKinDyn->getCenterOfMassJacobian(m_pimpl->robotState, m_pimpl->comJacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
    assert(ok);

    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);
    iDynTree::iDynTreeEigenMatrixMap comJacobianMap = iDynTree::toEigen(m_pimpl->comJacobianBuffer);

    iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

    jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset) = comJacobianMap.topLeftCorner<3, 3>();
    jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) = comJacobianMap.block<3, 3>(0, 3) * iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap);
    jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) = comJacobianMap.topRightCorner(3, m_pimpl->jointsPositionRange.size);

    jacobianMap.block<3,3>(0, m_pimpl->comPositionRange.offset).setIdentity();
    jacobianMap.block<3,3>(0, m_pimpl->comPositionRange.offset) *= -1;

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;
}

bool CoMPositionConstraint::constraintJacobianWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &jacobian)
{
    jacobian.zero();
    return true;
}

size_t CoMPositionConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t CoMPositionConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}
