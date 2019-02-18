/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Costs.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
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
        ok = stateVariables.addLabel(footName + "PositionPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
        ok = controlVariables.addLabel(footName + "VelocityControlPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
        ok = controlVariables.addLabel(footName + "ForceControlPoint" + std::to_string(i), 3);
        ASSERT_IS_TRUE(ok);
    }
}

void setVariables(VariablesLabeller& stateVariables, VariablesLabeller& controlVariables, size_t numberOfPoints, std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn) {

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

    ok = stateVariables.addLabel("JointsPosition", timelySharedKinDyn->model().getNrOfDOFs());
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("BaseLinearVelocity", 3);
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("BaseQuaternionDerivative", 4);
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("JointsVelocity", timelySharedKinDyn->model().getNrOfDOFs());
    ASSERT_IS_TRUE(ok);
}

void configureSharedKinDyn(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn) {
    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});
//    std::vector<std::string> vectorList({"r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch"});
    iDynTree::ModelLoader modelLoader;
    bool ok = modelLoader.loadModelFromFile(getAbsModelPath("iCubGenova04.urdf"));
    ASSERT_IS_TRUE(ok);
    ok = modelLoader.loadReducedModelFromFullModel(modelLoader.model(), vectorList);
    ASSERT_IS_TRUE(ok);
    assert(timelySharedKinDyn);
    ok = timelySharedKinDyn->loadRobotModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
//    ASSERT_IS_TRUE(sharedKinDyn->model().getNrOfDOFs() == 23);
    std::vector<double> timings(2);
    timings[0] = 0.0;
    timings[1] = 1.0;

    ok = timelySharedKinDyn->setTimings(timings);
    ASSERT_IS_TRUE(ok);

}

void configureCosts(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                    std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, std::shared_ptr<ExpressionsServer> expressionsServer,
                    const std::vector<iDynTree::Position>& leftPositions, const std::vector<iDynTree::Position>& rightPositions,
                    std::shared_ptr<StaticTorquesCost>& staticTorquesCost,
                    iDynTree::optimalcontrol::OptimalControlProblem& ocProblem) {
    bool ok = false;
    HyperbolicSecant forceActivation;
    forceActivation.setScaling(1.0);

    for (size_t i = 0; i < leftPositions.size(); ++i) {
        std::shared_ptr<ForceMeanCost> forceCost = std::make_shared<ForceMeanCost>(stateVariables, controlVariables, "Left", i);
        ok = ocProblem.addLagrangeTerm(1.0, forceCost);
        ASSERT_IS_TRUE(ok);

        std::shared_ptr<SwingCost> swingCost = std::make_shared<SwingCost>(stateVariables, controlVariables, "Left", i, 0.03);
        ok = ocProblem.addLagrangeTerm(1.0, swingCost);
        ASSERT_IS_TRUE(ok);

        std::shared_ptr<PhantomForcesCost> phantomForceCost = std::make_shared<PhantomForcesCost>(stateVariables, controlVariables, "Left",
                                                                                                  i, forceActivation);
        ok = ocProblem.addLagrangeTerm(1.0, phantomForceCost);
        ASSERT_IS_TRUE(ok);
    }
    for (size_t i = 0; i < rightPositions.size(); ++i) {
        std::shared_ptr<ForceMeanCost> forceCost = std::make_shared<ForceMeanCost>(stateVariables, controlVariables, "Right", i);
        ok = ocProblem.addLagrangeTerm(1.0, forceCost);
        ASSERT_IS_TRUE(ok);

        std::shared_ptr<SwingCost> swingCost = std::make_shared<SwingCost>(stateVariables, controlVariables, "Right", i, 0.03);
        ok = ocProblem.addLagrangeTerm(1.0, swingCost);
        ASSERT_IS_TRUE(ok);

        std::shared_ptr<PhantomForcesCost> phantomForceCost = std::make_shared<PhantomForcesCost>(stateVariables, controlVariables, "Right",
                                                                                                  i, forceActivation);
        ok = ocProblem.addLagrangeTerm(1.0, phantomForceCost);
        ASSERT_IS_TRUE(ok);
    }
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> comCost = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("CoMCost", stateVariables.getIndexRange("CoMPosition"),
                                                                                                                           stateVariables.size(),iDynTree::IndexRange::InvalidRange(),
                                                                                                                           controlVariables.size());
    ok = ocProblem.addLagrangeTerm(0.5, comCost);
    ASSERT_IS_TRUE(ok);
    std::shared_ptr<FrameOrientationCost> orientationCost = std::make_shared<FrameOrientationCost>(stateVariables, controlVariables, timelySharedKinDyn, expressionsServer, 22);
    ok = ocProblem.addLagrangeTerm(0.5, orientationCost);
    ASSERT_IS_TRUE(ok);

    iDynTree::FrameIndex leftFrame = timelySharedKinDyn->model().getFrameIndex("l_sole"),
            rightFrame = timelySharedKinDyn->model().getFrameIndex("r_sole");
    staticTorquesCost = std::make_shared<StaticTorquesCost>(stateVariables, controlVariables, timelySharedKinDyn, leftFrame, rightFrame, leftPositions, rightPositions);
    ok = ocProblem.addLagrangeTerm(0.5, staticTorquesCost);
    ASSERT_IS_TRUE(ok);

    std::shared_ptr<MeanPointPositionCost> meanPositionCost = std::make_shared<MeanPointPositionCost>(stateVariables, controlVariables);
    meanPositionCost->setTimeVaryingWeight(std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(15.0));
    ok = ocProblem.addLagrangeTerm(0.5, meanPositionCost);
    ASSERT_IS_TRUE(ok);

}

void checkCostsDerivative(double time, const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
                                double perturbation, iDynTree::optimalcontrol::OptimalControlProblem &ocProblem) {
    double originalCost, perturbedCost, firstOrderTaylor;

    iDynTree::VectorDynSize perturbedState, perturbedControl;
    iDynTree::VectorDynSize stateJacobian, controlJacobian;

    bool ok = false;

    ok = ocProblem.costsEvaluation(time, originalStateVector, originalControlVector, originalCost);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.costsFirstPartialDerivativeWRTState(time, originalStateVector, originalControlVector, stateJacobian);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.costsFirstPartialDerivativeWRTControl(time, originalStateVector, originalControlVector, controlJacobian);
    ASSERT_IS_TRUE(ok);

    for (unsigned int i = 0; i < originalStateVector.size(); ++i) {
        perturbedState = originalStateVector;
        perturbedState(i) = perturbedState(i) + perturbation;

        ok = ocProblem.costsEvaluation(time, perturbedState, originalControlVector, perturbedCost);
        ASSERT_IS_TRUE(ok);

        firstOrderTaylor = originalCost + iDynTree::toEigen(stateJacobian).transpose() * (iDynTree::toEigen(perturbedState) - iDynTree::toEigen(originalStateVector));
        ASSERT_EQUAL_DOUBLE_TOL(perturbedCost, firstOrderTaylor, perturbation/10);
    }

    for (unsigned int i = 0; i < originalControlVector.size(); ++i) {
        perturbedControl = originalControlVector;
        perturbedControl(i) = perturbedControl(i) + perturbation;

        ok = ocProblem.costsEvaluation(time, originalStateVector, perturbedControl, perturbedCost);
        ASSERT_IS_TRUE(ok);

        firstOrderTaylor = originalCost + iDynTree::toEigen(controlJacobian).transpose() * (iDynTree::toEigen(perturbedControl) - iDynTree::toEigen(originalControlVector));
        ASSERT_EQUAL_DOUBLE_TOL(perturbedCost, firstOrderTaylor, perturbation/10);
    }
}

void checkFrameOrientationDerivative(const iDynTree::Rotation& desiredRotation, const iDynTree::FrameIndex& testFrame, double perturbation, std::shared_ptr<SharedKinDynComputations> sharedKinDyn, const VariablesLabeller& stateVariables) {
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

    robotState.base_quaternion = baseQuaternion;
    robotState.base_position = basePosition;

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

        robotState.base_quaternion = baseQuaternion;
        robotState.base_position = basePosition;

        robotState.s = perturbedState(jointsPositionRange);


        frameTransform = sharedKinDyn->getWorldTransform(robotState, testFrame);
        quaternionErrorPerturbed = ErrorQuaternion(frameTransform.getRotation(), desiredRotation);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(quaternionError) +
                jacobianMap * (iDynTree::toEigen(perturbedState.values()) - iDynTree::toEigen(stateVariables.values()));
        ASSERT_EQUAL_VECTOR_TOL(quaternionErrorPerturbed, firstOrderTaylor, perturbation/10);
    }
}

void checkStaticForcesJacobian(double time, const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
        double perturbation, std::shared_ptr<StaticTorquesCost> staticTorques, unsigned int initialIndex) {
    iDynTree::VectorDynSize originalTorques, perturbedTorques, perturbedState, perturbedControl, firstOrderTaylor;
    iDynTree::MatrixDynSize stateJacobian;


    staticTorques->computeStaticTorques(time, originalStateVector, originalControlVector, originalTorques);

    perturbedTorques = originalTorques;
    firstOrderTaylor = originalTorques;

    staticTorques->computeStaticTorquesJacobian(time, originalStateVector, originalControlVector, stateJacobian);

//    std::cerr << "Jacobian: " << std::endl << stateJacobian.toString() << std::endl;


    for (unsigned int i = initialIndex; i < originalStateVector.size(); ++i) {
//        std::cerr << "State: " << i <<std::endl;
        perturbedState = originalStateVector;
        perturbedState(i) = perturbedState(i) + perturbation;

        staticTorques->computeStaticTorques(time, perturbedState, originalControlVector, perturbedTorques);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalTorques) +
                iDynTree::toEigen(stateJacobian) * (iDynTree::toEigen(perturbedState) - iDynTree::toEigen(originalStateVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedTorques, firstOrderTaylor, perturbation/10);
    }
}


int main() {
    VariablesLabeller stateVariables, controlVariables;
    std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn = std::make_shared<TimelySharedKinDynComputations>();
    iDynTree::optimalcontrol::OptimalControlProblem ocProblem;

    std::vector<iDynTree::Position> leftPositions, rightPositions;

    leftPositions.push_back(iDynTree::Position(0.125, -0.04, 0.0));
    leftPositions.push_back(iDynTree::Position(0.125,  0.04, 0.0));
    leftPositions.push_back(iDynTree::Position(-0.063,  0.04, 0.0));
    leftPositions.push_back(iDynTree::Position( 0.063, -0.04, 0.0));

    rightPositions.push_back(iDynTree::Position(0.125,  0.04, 0.0));
    rightPositions.push_back(iDynTree::Position(0.125, -0.04, 0.0));
    rightPositions.push_back(iDynTree::Position(-0.063, -0.04, 0.0));
    rightPositions.push_back(iDynTree::Position( 0.063,  0.04, 0.0));

    configureSharedKinDyn(timelySharedKinDyn);
    std::shared_ptr<ExpressionsServer> expressionsServer = std::make_shared<ExpressionsServer>(timelySharedKinDyn);

    setVariables(stateVariables, controlVariables, leftPositions.size(), timelySharedKinDyn);


    std::shared_ptr<StaticTorquesCost> staticTorquesPtr;

    configureCosts(stateVariables, controlVariables, timelySharedKinDyn, expressionsServer, leftPositions, rightPositions, staticTorquesPtr, ocProblem);

    iDynTree::VectorDynSize stateVector, controlVector;
    stateVector.resize(static_cast<unsigned int>(stateVariables.size()));
    iDynTree::getRandomVector(stateVector);
    controlVector.resize(static_cast<unsigned int>(controlVariables.size()));
    iDynTree::getRandomVector(controlVector);

    checkStaticForcesJacobian(0.0, stateVector, controlVector, 0.001, staticTorquesPtr, 0);

    checkCostsDerivative(0.0, stateVector, controlVector, 0.0001, ocProblem);
    checkCostsDerivative(1.0, stateVector, controlVector, 0.0001, ocProblem);


    stateVariables = stateVector;

    iDynTree::Rotation desiredRotation = iDynTree::getRandomRotation();//iDynTree::Rotation::Identity();

    checkFrameOrientationDerivative(desiredRotation, 0, 0.01, timelySharedKinDyn->get(0.0), stateVariables);
}
