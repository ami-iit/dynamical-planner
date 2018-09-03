/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#include <DynamicalPlannerPrivate/ContactPositionConsistencyConstraint.h>
#include <iDynTree/KinDynComputations.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/QuaternionUtils.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class ContactPositionConsistencyConstraint::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;

    std::string footName;
    iDynTree::FrameIndex footFrame;
    size_t contactIndex;
    iDynTree::Position positionInFoot;

    iDynTree::IndexRange positionPointRange, velocityPointRange;
    iDynTree::IndexRange jointsPositionRange, basePositionRange, baseQuaternionRange;

    iDynTree::VectorDynSize constraintValueBuffer;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::MatrixDynSize footJacobianBuffer, pointJacobianBuffer, stateJacobianBuffer, controlJacobianBuffer;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;
    iDynTree::MatrixFixSize<3, 6> footInternalTransformation, footTransformationBuffer;

    iDynTree::Rotation baseRotation;
    iDynTree::Position basePosition;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputation> sharedKinDyn;

    void getRanges() {

        velocityPointRange = stateVariables.getIndexRange(footName + "VelocityPoint" + std::to_string(contactIndex));
        assert(velocityPointRange.isValid());

        positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
        assert(positionPointRange.isValid());

        basePositionRange = stateVariables.getIndexRange("BasePosition");
        assert(basePositionRange.isValid());

        baseQuaternionRange = stateVariables.getIndexRange("BaseQuaternion");
        assert(baseQuaternionRange.isValid());

        jointsPositionRange = stateVariables.getIndexRange("JointsPosition");
        assert(jointsPositionRange.isValid());
    }

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



ContactPositionConsistencyConstraint::ContactPositionConsistencyConstraint(const VariablesLabeller& stateVariables,
                                                                           const VariablesLabeller& controlVariables,
                                                                           std::shared_ptr<SharedKinDynComputation> sharedKinDyn,
                                                                           iDynTree::FrameIndex footFrame, const std::string &footName,
                                                                           const iDynTree::Position &positionInFoot, size_t contactIndex)
    : iDynTree::optimalcontrol::Constraint(3, "ContactPositionConsistency" + footName + std::to_string(contactIndex))
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    assert(sharedKinDyn);
    assert(sharedKinDyn->isValid());

    m_pimpl->sharedKinDyn = sharedKinDyn;
    m_pimpl->footFrame = footFrame;
    assert(footFrame != iDynTree::FRAME_INVALID_INDEX);
    m_pimpl->footName = footName;
    m_pimpl->positionInFoot = positionInFoot;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->getRanges();

    m_pimpl->constraintValueBuffer.resize(3);
    m_pimpl->constraintValueBuffer.zero();
    m_pimpl->footJacobianBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->pointJacobianBuffer.resize(3, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateJacobianBuffer.resize(3, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();
    m_pimpl->controlJacobianBuffer.resize(3, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->robotState = sharedKinDyn->currentState();

    iDynTree::toEigen(m_pimpl->footInternalTransformation).leftCols<3>().setIdentity();
    iDynTree::toEigen(m_pimpl->footInternalTransformation).rightCols<3>() = iDynTree::skew(-iDynTree::toEigen(m_pimpl->positionInFoot));

    m_isLowerBounded = true;
    m_isUpperBounded = true;

}

ContactPositionConsistencyConstraint::~ContactPositionConsistencyConstraint()
{ }

bool ContactPositionConsistencyConstraint::evaluateConstraint(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;

    m_pimpl->updateRobotState();

    iDynTree::toEigen(m_pimpl->constraintValueBuffer) = iDynTree::toEigen(m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame) * m_pimpl->positionInFoot)
                                                            - iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->positionPointRange));

    constraint = m_pimpl->constraintValueBuffer;

    return true;
}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;

    m_pimpl->updateRobotState();

    bool ok = m_pimpl->sharedKinDyn->getFrameFreeFloatingJacobian(m_pimpl->robotState, m_pimpl->footFrame, m_pimpl->footJacobianBuffer, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
    assert(ok);

    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);

    Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::RowMajor> > footTransformationMap = iDynTree::toEigen(m_pimpl->footTransformationBuffer);

    footTransformationMap = iDynTree::toEigen(m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame).getRotation()) * iDynTree::toEigen(m_pimpl->footInternalTransformation);

    iDynTree::iDynTreeEigenMatrixMap footJacobianMap = iDynTree::toEigen(m_pimpl->footJacobianBuffer);
    iDynTree::iDynTreeEigenMatrixMap pointJacobianMap = iDynTree::toEigen(m_pimpl->pointJacobianBuffer);

    iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

    pointJacobianMap = footTransformationMap * footJacobianMap;

    jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset) = pointJacobianMap.leftCols<3>();
    jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) = pointJacobianMap.block<3, 3>(0, 3) * iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap);
    jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) = pointJacobianMap.rightCols(m_pimpl->jointsPositionRange.size);

    jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset).setIdentity();
    jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset) *= -1;

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;

}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ContactPositionConsistencyConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ContactPositionConsistencyConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}
