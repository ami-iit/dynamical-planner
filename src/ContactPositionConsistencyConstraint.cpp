/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#include <levi/levi.h>
#include <DynamicalPlannerPrivate/ContactPositionConsistencyConstraint.h>
#include <DynamicalPlannerPrivate/CheckEqualVector.h>
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

    iDynTree::IndexRange positionPointRange;
    iDynTree::IndexRange jointsPositionRange, basePositionRange, baseQuaternionRange;

    iDynTree::VectorDynSize constraintValueBuffer;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::Position pointPosition;
    iDynTree::MatrixDynSize footJacobianBuffer, pointJacobianBuffer, stateJacobianBuffer, controlJacobianBuffer;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;
    iDynTree::MatrixFixSize<3, 6> /*footInternalTransformation,*/ footTransformationBuffer;

    iDynTree::Rotation baseRotation;
    iDynTree::Position basePosition;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;

    bool updateDoneOnceConstraint = false;
    bool updateDoneOnceStateJacobian = false;
    double tolerance;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;

    void getRanges() {

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
        baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);

        robotState.world_T_base.setRotation(baseRotation);
        robotState.world_T_base.setPosition(basePosition);

        robotState.s = stateVariables(jointsPositionRange);
    }

    void updateVariables (){
        updateRobotState();
        iDynTree::toEigen(pointPosition) = iDynTree::toEigen(stateVariables(positionPointRange));
    }

    bool sameVariables(bool updateDoneOnce) {
        bool same = updateDoneOnce;
        same = same && VectorsAreEqual(pointPosition, stateVariables(positionPointRange), tolerance);
        same = same && VectorsAreEqual(basePosition, stateVariables(basePositionRange), tolerance);
        same = same && VectorsAreEqual(baseQuaternion, stateVariables(baseQuaternionRange), tolerance);
        same = same && VectorsAreEqual(robotState.s, stateVariables(jointsPositionRange), tolerance);

        return same;
    }

    void setSparsity() {
        stateSparsity.clear();
        controlSparsity.clear();

        iDynTree::IndexRange fullRange;
        fullRange.offset = 0;
        fullRange.size = 3;

        stateSparsity.addIdentityBlock(0, static_cast<size_t>(basePositionRange.offset), 3);
        stateSparsity.addDenseBlock(fullRange, baseQuaternionRange);
        stateSparsity.addDenseBlock(fullRange, jointsPositionRange);
        stateSparsity.addIdentityBlock(0, static_cast<size_t>(positionPointRange.offset), 3);

    }
};



ContactPositionConsistencyConstraint::ContactPositionConsistencyConstraint(const VariablesLabeller& stateVariables,
                                                                           const VariablesLabeller& controlVariables,
                                                                           std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                                                           iDynTree::FrameIndex footFrame, const std::string &footName,
                                                                           const iDynTree::Position &positionInFoot, size_t contactIndex)
    : iDynTree::optimalcontrol::Constraint(3, "ContactPositionConsistency" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());

    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;
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

    m_pimpl->tolerance = timelySharedKinDyn->getUpdateTolerance();

    iDynTree::toEigen(m_pimpl->footTransformationBuffer).leftCols<3>().setIdentity();


    m_isLowerBounded = true;
    m_isUpperBounded = true;
    m_upperBound.zero();
    m_lowerBound.zero();

    m_pimpl->setSparsity();

}

ContactPositionConsistencyConstraint::~ContactPositionConsistencyConstraint()
{ }

void ContactPositionConsistencyConstraint::setEqualityTolerance(double tolerance)
{
    assert(tolerance >= 0);

    iDynTree::toEigen(m_lowerBound).setConstant(-tolerance/2.0);
    iDynTree::toEigen(m_upperBound).setConstant(tolerance/2.0);
}

bool ContactPositionConsistencyConstraint::evaluateConstraint(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceConstraint))) {

//        m_pimpl->updateDoneOnceConstraint = true;
        m_pimpl->updateVariables();

        iDynTree::toEigen(m_pimpl->constraintValueBuffer) = iDynTree::toEigen(m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame) * m_pimpl->positionInFoot
                - m_pimpl->pointPosition);
    }

    constraint = m_pimpl->constraintValueBuffer;

    return true;
}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceStateJacobian))) {

//        m_pimpl->updateDoneOnceStateJacobian = true;
        m_pimpl->updateVariables();

        bool ok = m_pimpl->sharedKinDyn->getFrameFreeFloatingJacobian(m_pimpl->robotState, m_pimpl->footFrame, m_pimpl->footJacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
        assert(ok);

        iDynTree::iDynTreeEigenMatrixMap footJacobianMap = iDynTree::toEigen(m_pimpl->footJacobianBuffer);
        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);
        iDynTree::Transform footTransform = m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame);

//        Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::RowMajor> > footTransformationMap = iDynTree::toEigen(m_pimpl->footTransformationBuffer);

//        //    footTransformationMap = iDynTree::toEigen(m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame).getRotation()) * iDynTree::toEigen(m_pimpl->footInternalTransformation);
//        footTransformationMap.rightCols<3>() = iDynTree::skew(iDynTree::toEigen(footTransform.getPosition() - (footTransform * m_pimpl->positionInFoot)));

//        iDynTree::iDynTreeEigenMatrixMap pointJacobianMap = iDynTree::toEigen(m_pimpl->pointJacobianBuffer);

//        //    iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));
//        iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

//        pointJacobianMap = footTransformationMap * footJacobianMap;

//        jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset) = pointJacobianMap.leftCols<3>();

//        jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) = pointJacobianMap.block<3, 3>(0, 3) * iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap);

//        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) = pointJacobianMap.rightCols(m_pimpl->jointsPositionRange.size);

//        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset).setIdentity();
//        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset) *= -1;


        iDynTree::Transform footTransformInBase = m_pimpl->robotState.world_T_base.inverse() * footTransform;
        iDynTree::Position expectedPointPositionInBase = footTransformInBase * m_pimpl->positionInFoot;
        iDynTree::Vector4 footQuaternion = footTransform.getRotation().asQuaternion();

        jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset).setIdentity();

        jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) =
                iDynTree::toEigen(RotatedVectorQuaternionJacobian(expectedPointPositionInBase, m_pimpl->baseQuaternionNormalized)) *
                iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) =
                    footJacobianMap.topRightCorner(3, m_pimpl->jointsPositionRange.size) +
                    (iDynTree::toEigen(RotatedVectorQuaternionJacobian(m_pimpl->positionInFoot, footQuaternion)) *
                     iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivative(footQuaternion)) *
                     footJacobianMap.bottomRightCorner(3, m_pimpl->jointsPositionRange.size));


        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset).setIdentity();
        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset) *= -1;
    }

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

bool ContactPositionConsistencyConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}
