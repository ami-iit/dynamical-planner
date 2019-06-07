/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Costs/FrameAngularVelocityCost.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/CheckEqualVector.h>
#include <cassert>
#include <iostream>

using namespace DynamicalPlanner::Private;

class FrameAngularVelocityCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;
    iDynTree::FrameIndex desiredFrame;

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;
    iDynTree::IndexRange jointsPositionRange, basePositionRange, baseQuaternionRange, baseQuaternionDerivativeRange, jointsVelocityRange;
    iDynTree::Position basePosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::Rotation baseRotation, rotationError;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    levi::Mutable desiredRotation = levi::Mutable(3,3,"DesiredRotation");

    levi::Expression asExpression, quaternionDerivative, jointsDerivative, baseQuaternionVelocityDerivative, jointsVelocityDerivative,
        quaternionHessian, jointsHessian, quaternionJointsHessian, quaternionQuaternionDerivativeHessian,
        quaternionJointsVelocityHessian, jointsQuaternionDerivativeHessian, jointsJointsVelocityHessian;

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingRotation> desiredTrajectory;

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

    void constructExpressions(double rotationalPIDgain) {
        std::string desiredFrameName = timedSharedKinDyn->model().getFrameName(desiredFrame);

        levi::Expression frameAngularVelocity = expressionsServer->absoluteVelocity(desiredFrameName).block(3,0,3,1);
        levi::Expression frameRotation = expressionsServer->baseRotation() *
            expressionsServer->relativeRotation(timedSharedKinDyn->getFloatingBase(), desiredFrameName);

        levi::Expression rotationErrorExpression = frameRotation * desiredRotation.transpose();

        levi::Expression skewError = 0.5 * (rotationErrorExpression - rotationErrorExpression.transpose());

        levi::Expression velocityError = frameAngularVelocity + rotationalPIDgain * skewError.vee();

        asExpression = 0.5 * velocityError.transpose() * velocityError;

        quaternionDerivative = (velocityError.transpose() * velocityError.getColumnDerivative(0, expressionsServer->baseQuaternion())).transpose();
        jointsDerivative = (velocityError.transpose() * velocityError.getColumnDerivative(0, expressionsServer->jointsPosition())).transpose();
        baseQuaternionVelocityDerivative = (velocityError.transpose() * velocityError.getColumnDerivative(0, expressionsServer->baseQuaternionVelocity())).transpose();
        jointsVelocityDerivative = (velocityError.transpose() * velocityError.getColumnDerivative(0, expressionsServer->jointsVelocity())).transpose();

        quaternionHessian = quaternionDerivative.getColumnDerivative(0, expressionsServer->baseQuaternion());
        jointsHessian = jointsDerivative.getColumnDerivative(0, expressionsServer->jointsPosition());
        quaternionJointsHessian = quaternionDerivative.getColumnDerivative(0, expressionsServer->jointsPosition());
        quaternionQuaternionDerivativeHessian = quaternionDerivative.getColumnDerivative(0, expressionsServer->baseQuaternionVelocity());
        quaternionJointsVelocityHessian = quaternionDerivative.getColumnDerivative(0, expressionsServer->jointsVelocity());
        jointsQuaternionDerivativeHessian = jointsDerivative.getColumnDerivative(0, expressionsServer->baseQuaternionVelocity());
        jointsJointsVelocityHessian = jointsDerivative.getColumnDerivative(0, expressionsServer->jointsVelocity());
    }

};

FrameAngularVelocityCost::FrameAngularVelocityCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                   std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                                   std::shared_ptr<ExpressionsServer> expressionsServer,
                                                   const iDynTree::FrameIndex &desiredFrame, double rotationalPIDgain)
    : iDynTree::optimalcontrol::Cost ("FrameAngularVelocity")
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

    m_pimpl->baseQuaternionDerivativeRange = m_pimpl->controlVariables.getIndexRange("BaseQuaternionDerivative");
    assert(m_pimpl->baseQuaternionDerivativeRange.isValid());

    m_pimpl->jointsVelocityRange = m_pimpl->controlVariables.getIndexRange("JointsVelocity");
    assert(m_pimpl->jointsVelocityRange.isValid());

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();

    m_pimpl->expressionsServer = expressionsServer;

    m_pimpl->desiredTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantRotation>(iDynTree::Rotation::Identity());

    m_pimpl->constructExpressions(rotationalPIDgain);

    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->baseQuaternionRange, m_pimpl->baseQuaternionRange);
    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->baseQuaternionRange, m_pimpl->jointsPositionRange);
    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->jointsPositionRange, m_pimpl->baseQuaternionRange);
    m_pimpl->stateHessianSparsity.addDenseBlock(m_pimpl->jointsPositionRange, m_pimpl->jointsPositionRange);

    m_pimpl->mixedHessianSparsity.addDenseBlock(m_pimpl->baseQuaternionRange, m_pimpl->baseQuaternionDerivativeRange);
    m_pimpl->mixedHessianSparsity.addDenseBlock(m_pimpl->baseQuaternionRange, m_pimpl->jointsVelocityRange);
    m_pimpl->mixedHessianSparsity.addDenseBlock(m_pimpl->jointsPositionRange, m_pimpl->baseQuaternionDerivativeRange);
    m_pimpl->mixedHessianSparsity.addDenseBlock(m_pimpl->jointsPositionRange, m_pimpl->jointsVelocityRange);

    m_pimpl->controlHessianSparsity.clear();

}

FrameAngularVelocityCost::~FrameAngularVelocityCost()
{
    m_pimpl->asExpression.clearDerivativesCache();
    m_pimpl->quaternionDerivative.clearDerivativesCache();
    m_pimpl->jointsDerivative.clearDerivativesCache();
    m_pimpl->baseQuaternionVelocityDerivative.clearDerivativesCache();
    m_pimpl->jointsVelocityDerivative.clearDerivativesCache();
    m_pimpl->quaternionHessian.clearDerivativesCache();
    m_pimpl->quaternionJointsHessian.clearDerivativesCache();
    m_pimpl->jointsHessian.clearDerivativesCache();
    m_pimpl->quaternionQuaternionDerivativeHessian.clearDerivativesCache();
    m_pimpl->quaternionJointsVelocityHessian.clearDerivativesCache();
    m_pimpl->jointsQuaternionDerivativeHessian.clearDerivativesCache();
    m_pimpl->jointsJointsVelocityHessian.clearDerivativesCache();
}

bool FrameAngularVelocityCost::setDesiredRotationTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingRotation> desiredRotationTrajectory)
{
    if (!desiredRotationTrajectory)
        return false;
    m_pimpl->desiredTrajectory = desiredRotationTrajectory;
    return true;
}

bool FrameAngularVelocityCost::costEvaluation(double time, const iDynTree::VectorDynSize &state,
                                          const iDynTree::VectorDynSize &control, double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time);

    bool isValid = false;
    const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid rotation at time " << time
                  << "." << std::endl;
        return false;
    }

    m_pimpl->desiredRotation = iDynTree::toEigen(desiredRotation);

    costValue = m_pimpl->asExpression.evaluate()(0);

    return true;
}

bool FrameAngularVelocityCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                                  const iDynTree::VectorDynSize &control,
                                                                  iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time);

    bool isValid = false;
    const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid rotation at time " << time
                  << "." << std::endl;
        return false;
    }

    m_pimpl->desiredRotation = iDynTree::toEigen(desiredRotation);

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);

    gradientMap.segment<4>(m_pimpl->baseQuaternionRange.offset) = m_pimpl->quaternionDerivative.evaluate();
    gradientMap.segment(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.size) =
        m_pimpl->jointsDerivative.evaluate();

    partialDerivative = m_pimpl->stateGradientBuffer;

    return true;
}

bool FrameAngularVelocityCost::costFirstPartialDerivativeWRTControl(double time, const iDynTree::VectorDynSize &state,
                                                                    const iDynTree::VectorDynSize &control,
                                                                    iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time);

    bool isValid = false;
    const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid rotation at time " << time
                  << "." << std::endl;
        return false;
    }

    m_pimpl->desiredRotation = iDynTree::toEigen(desiredRotation);

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->controlGradientBuffer);

    gradientMap.segment<4>(m_pimpl->baseQuaternionDerivativeRange.offset) = m_pimpl->baseQuaternionVelocityDerivative.evaluate();
    gradientMap.segment(m_pimpl->jointsVelocityRange.offset, m_pimpl->jointsVelocityRange.size) =
        m_pimpl->jointsVelocityDerivative.evaluate();

    partialDerivative = m_pimpl->controlGradientBuffer;

    return true;
}

bool FrameAngularVelocityCost::costSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                                   const iDynTree::VectorDynSize &control,
                                                                   iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time);

    bool isValid = false;
    const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid rotation at time " << time
                  << "." << std::endl;
        return false;
    }

    m_pimpl->desiredRotation = iDynTree::toEigen(desiredRotation);

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

bool FrameAngularVelocityCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                 const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FrameAngularVelocityCost::costSecondPartialDerivativeWRTStateControl(double time, const iDynTree::VectorDynSize &state,
                                                                          const iDynTree::VectorDynSize &control,
                                                                          iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time);

    bool isValid = false;
    const iDynTree::Rotation& desiredRotation = m_pimpl->desiredTrajectory->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][FrameOrientationCost::costEvaluation] Unable to retrieve a valid rotation at time " << time
                  << "." << std::endl;
        return false;
    }

    m_pimpl->desiredRotation = iDynTree::toEigen(desiredRotation);

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(partialDerivative);

    hessianMap.block<4,4>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseQuaternionDerivativeRange.offset) =
        m_pimpl->quaternionQuaternionDerivativeHessian.evaluate();

    hessianMap.block(m_pimpl->baseQuaternionRange.offset, m_pimpl->jointsVelocityRange.offset, m_pimpl->baseQuaternionRange.size, m_pimpl->jointsVelocityRange.size) = m_pimpl->quaternionJointsVelocityHessian.evaluate();

    hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsVelocityRange.offset,
                     m_pimpl->jointsPositionRange.size, m_pimpl->jointsVelocityRange.size) =
        m_pimpl->jointsJointsVelocityHessian.evaluate();


    hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionDerivativeRange.offset, m_pimpl->jointsPositionRange.size, m_pimpl->baseQuaternionDerivativeRange.size) = m_pimpl->jointsQuaternionDerivativeHessian.evaluate();

    return true;
}

bool FrameAngularVelocityCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool FrameAngularVelocityCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool FrameAngularVelocityCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
