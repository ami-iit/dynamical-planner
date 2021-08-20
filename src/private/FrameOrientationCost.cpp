/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Costs/FrameOrientationCost.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/CheckEqualVector.h>
#include <cassert>
#include <iostream>

using namespace DynamicalPlanner::Private;

class FrameOrientationCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;
    iDynTree::FrameIndex desiredFrame;

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer, reducedGradientBuffer;
    iDynTree::IndexRange jointsPositionRange, basePositionRange, baseQuaternionRange;
    iDynTree::MatrixDynSize frameJacobianBuffer, frameJacobianTimesMap, quaternionErrorPartialDerivative;
    iDynTree::Position basePosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized, identityQuaternion, quaternionError, quaternionDifference;
    iDynTree::Rotation baseRotation, rotationError;
    iDynTree::Transform frameTransform;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    levi::Expression quaternionErrorExpression,
        asExpression, quaternionDerivative, jointsDerivative,
        quaternionHessian, jointsHessian, quaternionJointsHessian;
    levi::Variable desiredQuaternion, timeWeightVariable;

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingRotation> desiredTrajectory;
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> timeVaryingWeight;

    bool updateDoneOnceCost = false;
    bool updateDoneOnceStateJacobian = false;
    double lastUpdateTimeCost = -1;
    double lastUpdateTimeStateJacobian = -1;
    double costBuffer;
    double tolerance;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

    void updateRobotState() {

        robotState = sharedKinDyn->currentState();

        iDynTree::toEigen(basePosition) = iDynTree::toEigen(stateVariables(basePositionRange));
        baseQuaternion = stateVariables(baseQuaternionRange);
        baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);

        robotState.base_quaternion = baseQuaternion;
        robotState.base_position = basePosition;

        robotState.s = stateVariables(jointsPositionRange);
        sharedKinDyn->updateRobotState(robotState);

    }

    bool sameVariables(bool updateDoneOnce, double newTime, double lastUpdateTime) {
        bool same = updateDoneOnce;
        same = same && iDynTree::checkDoublesAreEqual(lastUpdateTime, newTime, tolerance);
        same = same && VectorsAreEqual(basePosition, stateVariables(basePositionRange), tolerance);
        same = same && VectorsAreEqual(baseQuaternion, stateVariables(baseQuaternionRange), tolerance);
        same = same && VectorsAreEqual(robotState.s, stateVariables(jointsPositionRange), tolerance);

        return same;
    }

};

FrameOrientationCost::FrameOrientationCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                           std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                           std::shared_ptr<ExpressionsServer> expressionsServer, const iDynTree::FrameIndex &desiredFrame)
    : iDynTree::optimalcontrol::Cost ("FrameOrientation")
    , m_pimpl(std::make_unique<Implementation>())
{
    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());
    assert(desiredFrame != iDynTree::FRAME_INVALID_INDEX);
    assert(timelySharedKinDyn->model().isValidFrameIndex(desiredFrame));

    m_pimpl->desiredFrame = desiredFrame;

    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->frameJacobianBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->frameJacobianTimesMap.resize(3, 7 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->frameJacobianTimesMap.zero();
    m_pimpl->quaternionErrorPartialDerivative.resize(4, 7 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();
    m_pimpl->reducedGradientBuffer.resize(7 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->reducedGradientBuffer.zero();

    m_pimpl->desiredTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantRotation>(iDynTree::Rotation::Identity());
    m_pimpl->timeVaryingWeight = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(1.0);
    m_pimpl->identityQuaternion = iDynTree::Rotation::Identity().asQuaternion();

    m_pimpl->expressionsServer = expressionsServer;

    std::string desiredFrameName = timelySharedKinDyn->model().getFrameName(desiredFrame);
    m_pimpl->desiredQuaternion = levi::Variable(4, "quat_desired");
    m_pimpl->timeWeightVariable = levi::Variable(1, "timeVaryingWeight");
    m_pimpl->quaternionErrorExpression = expressionsServer->quaternionError(desiredFrameName, m_pimpl->desiredQuaternion);
    levi::Constant identityQuat_expr(iDynTree::toEigen(m_pimpl->identityQuaternion),"quaternionIdentity");

    levi::Expression quaternionDifference = m_pimpl->quaternionErrorExpression - identityQuat_expr;
    m_pimpl->asExpression = 0.5 * m_pimpl->timeWeightVariable *  quaternionDifference.transpose() * quaternionDifference;

    m_pimpl->quaternionDerivative = m_pimpl->timeWeightVariable * (quaternionDifference.transpose() * quaternionDifference.getColumnDerivative(0, m_pimpl->expressionsServer->baseQuaternion())).transpose();
    m_pimpl->jointsDerivative = m_pimpl->timeWeightVariable * (quaternionDifference.transpose() * quaternionDifference.getColumnDerivative(0, m_pimpl->expressionsServer->jointsPosition())).transpose();

    m_pimpl->quaternionHessian = m_pimpl->quaternionDerivative.getColumnDerivative(0, m_pimpl->expressionsServer->baseQuaternion());
    m_pimpl->jointsHessian = m_pimpl->jointsDerivative.getColumnDerivative(0, m_pimpl->expressionsServer->jointsPosition());
    m_pimpl->quaternionJointsHessian = m_pimpl->quaternionDerivative.getColumnDerivative(0, m_pimpl->expressionsServer->jointsPosition());

    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->baseQuaternionRange, m_pimpl->baseQuaternionRange);
    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->baseQuaternionRange, m_pimpl->jointsPositionRange);
    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->jointsPositionRange, m_pimpl->baseQuaternionRange);
    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->jointsPositionRange, m_pimpl->jointsPositionRange);

    m_pimpl->mixedHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();

}

FrameOrientationCost::~FrameOrientationCost()
{
    m_pimpl->quaternionDerivative.clearDerivativesCache();
    m_pimpl->jointsDerivative.clearDerivativesCache();
    m_pimpl->quaternionHessian.clearDerivativesCache();
    m_pimpl->quaternionJointsHessian.clearDerivativesCache();
    m_pimpl->jointsHessian.clearDerivativesCache();
}

void FrameOrientationCost::setDesiredRotation(const iDynTree::Rotation &desiredRotation)
{
    m_pimpl->desiredTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantRotation>(desiredRotation);
}

bool FrameOrientationCost::setDesiredRotationTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingRotation> desiredRotationTrajectory)
{
    if (!desiredRotationTrajectory)
        return false;
    m_pimpl->desiredTrajectory = desiredRotationTrajectory;
    return true;
}

void FrameOrientationCost::setTimeVaryingWeight(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> timeVaryingWeight)
{
    assert(timeVaryingWeight);
    m_pimpl->timeVaryingWeight = timeVaryingWeight;
}

bool FrameOrientationCost::costEvaluation(double time, const iDynTree::VectorDynSize &state,
                                          const iDynTree::VectorDynSize &control, double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceCost, time, m_pimpl->lastUpdateTimeCost))) {

//        m_pimpl->updateDoneOnceCost = true;
        m_pimpl->lastUpdateTimeCost = time;
        m_pimpl->updateRobotState();

        bool isValid = false;
        const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

        if (!isValid) {
            std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid rotation at time " << time
                      << "." << std::endl;
            return false;
        }

        const double& timeWeight = m_pimpl->timeVaryingWeight->get(time, isValid);

        if (!isValid) {
            std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid timeVaryingWeight at time " << time
                      << "." << std::endl;
            return false;
        }

        m_pimpl->frameTransform = m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->desiredFrame);

        m_pimpl->quaternionError = ErrorQuaternion(m_pimpl->frameTransform.getRotation(), desiredRotation);

        iDynTree::toEigen(m_pimpl->quaternionDifference) = iDynTree::toEigen(m_pimpl->quaternionError) -
                iDynTree::toEigen(m_pimpl->identityQuaternion);
        m_pimpl->costBuffer = 0.5 * timeWeight * QuaternionSquaredNorm(m_pimpl->quaternionDifference);
    }

    costValue = m_pimpl->costBuffer;

    return true;
}

bool FrameOrientationCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                              const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceStateJacobian, time, m_pimpl->lastUpdateTimeStateJacobian))) {

//        m_pimpl->updateDoneOnceStateJacobian = true;
        m_pimpl->lastUpdateTimeStateJacobian = time;
        m_pimpl->updateRobotState();

        bool isValid = false;
        const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

        if (!isValid) {
            std::cerr << "[ERROR][FrameOrientationCost::costFirstPartialDerivativeWRTState] Unable to retrieve a valid rotation at time " << time
                      << "." << std::endl;
            return false;
        }

        const double& timeWeight = m_pimpl->timeVaryingWeight->get(time, isValid);

        if (!isValid) {
            std::cerr << "[ERROR][FrameOrientationCost::costFirstPartialDerivativeWRTState] Unable to retrieve a valid timeVaryingWeight at time " << time
                      << "." << std::endl;
            return false;
        }

        m_pimpl->frameTransform = m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->desiredFrame);

        m_pimpl->quaternionError = ErrorQuaternion(m_pimpl->frameTransform.getRotation(), desiredRotation);

        iDynTree::toEigen(m_pimpl->quaternionDifference) = iDynTree::toEigen(m_pimpl->quaternionError) -
                iDynTree::toEigen(m_pimpl->identityQuaternion);

        bool ok = m_pimpl->sharedKinDyn->getFrameFreeFloatingJacobian(m_pimpl->robotState, m_pimpl->desiredFrame, m_pimpl->frameJacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

        assert(ok);

        iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) =
                iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) *
                iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

        iDynTree::iDynTreeEigenMatrixMap originalJacobian = iDynTree::toEigen(m_pimpl->frameJacobianBuffer);
        iDynTree::iDynTreeEigenMatrixMap jacobianModifiedMap = iDynTree::toEigen(m_pimpl->frameJacobianTimesMap);

        jacobianModifiedMap.leftCols<3>() = originalJacobian.block<3,3>(3,0);
        jacobianModifiedMap.block<3,4>(0, 3) = originalJacobian.block<3,3>(3,3) * iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap);
        jacobianModifiedMap.rightCols(m_pimpl->jointsPositionRange.size) = originalJacobian.block(3, 6, 3, m_pimpl->jointsPositionRange.size);

        iDynTree::iDynTreeEigenMatrixMap quaternionErrorDerivativeMap = iDynTree::toEigen(m_pimpl->quaternionErrorPartialDerivative);

        quaternionErrorDerivativeMap = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivative(m_pimpl->quaternionError)) *
                iDynTree::toEigen(desiredRotation).transpose() * jacobianModifiedMap;


        iDynTree::iDynTreeEigenVector reducedGradientMap = iDynTree::toEigen(m_pimpl->reducedGradientBuffer);

        reducedGradientMap = timeWeight * iDynTree::toEigen(m_pimpl->quaternionErrorPartialDerivative).transpose() * iDynTree::toEigen(m_pimpl->quaternionDifference);

        iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);

//        gradientMap.segment<3>(m_pimpl->basePositionRange.offset) = reducedGradientMap.segment<3>(0);
        gradientMap.segment<4>(m_pimpl->baseQuaternionRange.offset) = reducedGradientMap.segment<4>(3);
        gradientMap.segment(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.size) =
                reducedGradientMap.segment(7, m_pimpl->jointsPositionRange.size);

    }
    partialDerivative = m_pimpl->stateGradientBuffer;

    return true;
}

bool FrameOrientationCost::costFirstPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                const iDynTree::VectorDynSize &/*control*/,
                                                                iDynTree::VectorDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->controlGradientBuffer;

    return true;
}

bool FrameOrientationCost::costSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time);

    bool isValid = false;
    const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costSecondPartialDerivativeWRTState] Unable to retrieve a valid rotation at time " << time
                  << "." << std::endl;
        return false;
    }

    const double& timeWeight = m_pimpl->timeVaryingWeight->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costSecondPartialDerivativeWRTState] Unable to retrieve a valid timeVaryingWeight at time " << time
                  << "." << std::endl;
        return false;
    }

    m_pimpl->desiredQuaternion = iDynTree::toEigen(desiredRotation.asQuaternion());
    m_pimpl->timeWeightVariable = timeWeight;

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(partialDerivative);

    hessianMap.block<4,4>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseQuaternionRange.offset) =
        m_pimpl->quaternionHessian.evaluate();

    hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.size, m_pimpl->jointsPositionRange.size) =
        m_pimpl->jointsHessian.evaluate();


    const Eigen::MatrixXd& mixedHessian = m_pimpl->quaternionJointsHessian.evaluate();

    hessianMap.block(m_pimpl->baseQuaternionRange.offset, m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionRange.size, m_pimpl->jointsPositionRange.size) = mixedHessian;

    hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionRange.offset, m_pimpl->jointsPositionRange.size, m_pimpl->baseQuaternionRange.size) = mixedHessian.transpose();

    return true;
}

bool FrameOrientationCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                 const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FrameOrientationCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                      const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FrameOrientationCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool FrameOrientationCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool FrameOrientationCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
