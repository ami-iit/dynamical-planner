/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Costs.h>
#include <DynamicalPlannerPrivate/QuaternionUtils.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/OptimalControlProblem.h>
#include <URDFdir.h>
#include <memory>
#include <string>
#include <cassert>

using namespace DynamicalPlanner::Private;

void setFootVariables(VariablesLabeller& stateVariables, VariablesLabeller& controlVariables, const std::string& footName, size_t numberOfPoints) {
    bool ok = false;
    for (size_t i = 0; i < numberOfPoints; ++i) {
        ok = stateVariables.addLabel(footName + "ForcePoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
        ok = stateVariables.addLabel(footName + "VelocityPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
        ok = stateVariables.addLabel(footName + "PositionPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
        ok = controlVariables.addLabel(footName + "VelocityControlPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
        ok = controlVariables.addLabel(footName + "ForceControlPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
    }
}

void setVariables(VariablesLabeller& stateVariables, VariablesLabeller& controlVariables, size_t numberOfPoints) {

    setFootVariables(stateVariables, controlVariables, "Left", numberOfPoints);
    setFootVariables(stateVariables, controlVariables, "Right", numberOfPoints);

    bool ok;
    ok = stateVariables.addLabel("Momentum", 6);
    ASSERT_IS_TRUE(ok);

    ok = stateVariables.addLabel("CoMPosition", 3);
    ASSERT_IS_TRUE(ok);

    ok = stateVariables.addLabel("BasePosition", 3);
    ASSERT_IS_TRUE(ok);

    ok = stateVariables.addLabel("BaseQuaternion", 4);
    ASSERT_IS_TRUE(ok);

    ok = stateVariables.addLabel("JointsPosition", 23);
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("BaseVelocity", 6);
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("JointsVelocity", 23);
    ASSERT_IS_TRUE(ok);
}

void configureSharedKinDyn(std::shared_ptr<SharedKinDynComputation> sharedKinDyn) {
    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});
    iDynTree::ModelLoader modelLoader;
    bool ok = modelLoader.loadModelFromFile(getAbsModelPath("iCubGenova04.urdf"));
    ASSERT_IS_TRUE(ok);
    ok = modelLoader.loadReducedModelFromFullModel(modelLoader.model(), vectorList);
    ASSERT_IS_TRUE(ok);
    assert(sharedKinDyn);
    ok = sharedKinDyn->loadRobotModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
    ASSERT_IS_TRUE(sharedKinDyn->model().getNrOfDOFs() == 23);
}

void configureCosts(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                    std::shared_ptr<SharedKinDynComputation> sharedKinDyn, size_t numberOfPoints,
                    iDynTree::optimalcontrol::OptimalControlProblem& ocProblem) {
    bool ok = false;
    for (size_t i = 0; i < numberOfPoints; ++i) {
        std::shared_ptr<ForceMeanCost> forceCost = std::make_shared<ForceMeanCost>(stateVariables, controlVariables, "Left", i);
        ok = ocProblem.addLagrangeTerm(1.0, forceCost);
    }
    for (size_t i = 0; i < numberOfPoints; ++i) {
        std::shared_ptr<ForceMeanCost> forceCost = std::make_shared<ForceMeanCost>(stateVariables, controlVariables, "Right", i);
        ok = ocProblem.addLagrangeTerm(1.0, forceCost);
    }
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> comCost = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("CoMCost", stateVariables.getIndexRange("CoMPosition"),
                                                                                                                           stateVariables.size(),iDynTree::IndexRange::InvalidRange(),
                                                                                                                           controlVariables.size());
    ok = ocProblem.addLagrangeTerm(1.0, comCost);
    std::shared_ptr<FrameOrientationCost> orientationCost = std::make_shared<FrameOrientationCost>(stateVariables, controlVariables, sharedKinDyn, 22);
    ok = ocProblem.addLagrangeTerm(1.0, orientationCost);
}

void checkCostsDerivative(const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
                                double perturbation, iDynTree::optimalcontrol::OptimalControlProblem &ocProblem) {

    double originalCost, perturbedCost, firstOrderTaylor;

    iDynTree::VectorDynSize perturbedState, perturbedControl;
    iDynTree::VectorDynSize stateJacobian, controlJacobian;

    bool ok = false;

    ok = ocProblem.costsEvaluation(0.0, originalStateVector, originalControlVector, originalCost);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.costsFirstPartialDerivativeWRTState(0.0, originalStateVector, originalControlVector, stateJacobian);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.costsFirstPartialDerivativeWRTControl(0.0, originalStateVector, originalControlVector, controlJacobian);
    ASSERT_IS_TRUE(ok);

    for (unsigned int i = 0; i < originalStateVector.size(); ++i) {
        perturbedState = originalStateVector;
        perturbedState(i) = perturbedState(i) + perturbation;

        ok = ocProblem.costsEvaluation(0.0, perturbedState, originalControlVector, perturbedCost);
        ASSERT_IS_TRUE(ok);

        firstOrderTaylor = originalCost + iDynTree::toEigen(stateJacobian).transpose() * (iDynTree::toEigen(perturbedState) - iDynTree::toEigen(originalStateVector));
        ASSERT_EQUAL_DOUBLE_TOL(perturbedCost, firstOrderTaylor, perturbation/10);
    }

    for (unsigned int i = 0; i < originalControlVector.size(); ++i) {
        perturbedControl = originalControlVector;
        perturbedControl(i) = perturbedControl(i) + perturbation;

        ok = ocProblem.costsEvaluation(0.0, originalStateVector, perturbedControl, perturbedCost);
        ASSERT_IS_TRUE(ok);

        firstOrderTaylor = originalCost + iDynTree::toEigen(controlJacobian).transpose() * (iDynTree::toEigen(perturbedControl) - iDynTree::toEigen(originalControlVector));
        ASSERT_EQUAL_DOUBLE_TOL(perturbedCost, firstOrderTaylor, perturbation/10);
    }
}

void checkFrameOrientationDerivative(const iDynTree::Rotation& desiredRotation, const iDynTree::FrameIndex& testFrame, double perturbation, std::shared_ptr<SharedKinDynComputation> sharedKinDyn, const VariablesLabeller& stateVariables) {
    iDynTree::VectorDynSize perturbedControl, firstOrderTaylor(4);
    VariablesLabeller perturbedState = stateVariables;
    iDynTree::MatrixDynSize stateJacobian;

    iDynTree::IndexRange jointsPositionRange, basePositionRange, baseQuaternionRange;
    iDynTree::MatrixDynSize frameJacobianBuffer, frameJacobianTimesMap, quaternionErrorPartialDerivative, quaternionErrorPartialDerivativeTotal;
    iDynTree::Position basePosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized, quaternionError, quaternionErrorPerturbed, originalQuaternion;
    iDynTree::Rotation baseRotation, rotationError;
    iDynTree::Transform frameTransform;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;

    basePositionRange = stateVariables.getIndexRange("BasePosition");
    assert(basePositionRange.isValid());

    baseQuaternionRange = stateVariables.getIndexRange("BaseQuaternion");
    assert(baseQuaternionRange.isValid());

    jointsPositionRange = stateVariables.getIndexRange("JointsPosition");
    assert(jointsPositionRange.isValid());

    frameJacobianBuffer.resize(6, 6 + static_cast<unsigned int>(jointsPositionRange.size));
    frameJacobianTimesMap.resize(3, 7 + static_cast<unsigned int>(jointsPositionRange.size));
    frameJacobianTimesMap.zero();
    quaternionErrorPartialDerivative.resize(4, 7 + static_cast<unsigned int>(jointsPositionRange.size));
    quaternionErrorPartialDerivativeTotal.resize(4, static_cast<unsigned int>(stateVariables.size()));
    quaternionErrorPartialDerivativeTotal.zero();

    RobotState robotState;

    robotState = sharedKinDyn->currentState();

    iDynTree::toEigen(basePosition) = iDynTree::toEigen(stateVariables(basePositionRange));
    baseQuaternion = stateVariables(baseQuaternionRange);
    baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
    assert(QuaternionBoundsRespected(baseQuaternionNormalized));
    baseRotation.fromQuaternion(baseQuaternionNormalized);

    robotState.world_T_base.setRotation(baseRotation);
    robotState.world_T_base.setPosition(basePosition);

    robotState.s = stateVariables(jointsPositionRange);


    frameTransform = sharedKinDyn->getWorldTransform(robotState, testFrame);
    quaternionError = ErrorQuaternion(frameTransform.getRotation(), desiredRotation);

    bool ok = sharedKinDyn->getFrameFreeFloatingJacobian(robotState, testFrame, frameJacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    ASSERT_IS_TRUE(ok);

    iDynTree::toEigen(quaternionErrorPartialDerivative).leftCols<3>() = (iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivative(quaternionError)) * iDynTree::toEigen(desiredRotation.inverse()) *
                                                                         iDynTree::toEigen(frameJacobianBuffer).bottomRows<3>()).leftCols<3>();


    iDynTree::toEigen(quaternionErrorPartialDerivative).block<4,4>(0,3) = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivative(quaternionError)) * iDynTree::toEigen(desiredRotation.inverse()) * iDynTree::toEigen(frameJacobianBuffer).block<3,3>(3,3) *
            iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(baseQuaternion));


    iDynTree::toEigen(quaternionErrorPartialDerivative).rightCols(jointsPositionRange.size) = (iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivative(quaternionError)) * iDynTree::toEigen(desiredRotation.inverse()) *
                                                                                               iDynTree::toEigen(frameJacobianBuffer).bottomRows<3>()).rightCols(jointsPositionRange.size);

    // Left Trivialization case
//    bool ok = sharedKinDyn->getFrameFreeFloatingJacobian(robotState, testFrame, frameJacobianBuffer, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

//    ASSERT_IS_TRUE(ok);

//    std::cerr << "Frame Jacobian: " <<std::endl << frameJacobianBuffer.toString() << std::endl;

//    iDynTree::toEigen(quaternionErrorPartialDerivative).leftCols<3>() = (iDynTree::toEigen(QuaternionLeftTrivializedDerivative(quaternionError)) * iDynTree::toEigen(frameJacobianBuffer).bottomRows<3>()).leftCols<3>();


//    iDynTree::toEigen(quaternionErrorPartialDerivative).block<4,4>(0,3) = iDynTree::toEigen(QuaternionLeftTrivializedDerivative(quaternionError)) * iDynTree::toEigen(frameJacobianBuffer).block<3,3>(3,3) *
//            iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(baseQuaternion));


//    iDynTree::toEigen(quaternionErrorPartialDerivative).rightCols(jointsPositionRange.size) = (iDynTree::toEigen(QuaternionLeftTrivializedDerivative(quaternionError)) * iDynTree::toEigen(frameJacobianBuffer).bottomRows<3>()).rightCols(jointsPositionRange.size);


    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(quaternionErrorPartialDerivativeTotal);

    jacobianMap.block<4,3>(0, basePositionRange.offset) = iDynTree::toEigen(quaternionErrorPartialDerivative).block<4,3>(0,0);
    jacobianMap.block<4,4>(0, baseQuaternionRange.offset) = iDynTree::toEigen(quaternionErrorPartialDerivative).block<4,4>(0,3);
    jacobianMap.block(0, jointsPositionRange.offset, 4, jointsPositionRange.size) = iDynTree::toEigen(quaternionErrorPartialDerivative).block(0, 7, 4, jointsPositionRange.size);


    for (unsigned int i = 0; i < stateVariables.size(); ++i) {
        perturbedState = stateVariables;
        perturbedState(i) = perturbedState(i) + perturbation;

        robotState = sharedKinDyn->currentState();

        iDynTree::toEigen(basePosition) = iDynTree::toEigen(perturbedState(basePositionRange));
        baseQuaternion = perturbedState(baseQuaternionRange);
        baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);

        robotState.world_T_base.setRotation(baseRotation);
        robotState.world_T_base.setPosition(basePosition);

        robotState.s = perturbedState(jointsPositionRange);


        frameTransform = sharedKinDyn->getWorldTransform(robotState, testFrame);
        quaternionErrorPerturbed = ErrorQuaternion(frameTransform.getRotation(), desiredRotation);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(quaternionError) +
                jacobianMap * (iDynTree::toEigen(perturbedState.values()) - iDynTree::toEigen(stateVariables.values()));
        ASSERT_EQUAL_VECTOR_TOL(quaternionErrorPerturbed, firstOrderTaylor, perturbation/10);
    }


}


int main() {
    VariablesLabeller stateVariables, controlVariables;
    std::shared_ptr<SharedKinDynComputation> sharedKinDyn = std::make_shared<SharedKinDynComputation>();
    iDynTree::optimalcontrol::OptimalControlProblem ocProblem;

    setVariables(stateVariables, controlVariables, 4);
    configureSharedKinDyn(sharedKinDyn);
    configureCosts(stateVariables, controlVariables, sharedKinDyn, 4, ocProblem);

    iDynTree::VectorDynSize stateVector, controlVector;
    stateVector.resize(static_cast<unsigned int>(stateVariables.size()));
    iDynTree::getRandomVector(stateVector);
    controlVector.resize(static_cast<unsigned int>(controlVariables.size()));
    iDynTree::getRandomVector(controlVector);

    checkCostsDerivative(stateVector, controlVector, 0.001, ocProblem);

    stateVariables = stateVector;

    iDynTree::Rotation desiredRotation = iDynTree::getRandomRotation();//iDynTree::Rotation::Identity();

    checkFrameOrientationDerivative(desiredRotation, 0, 0.01, sharedKinDyn, stateVariables);
}
