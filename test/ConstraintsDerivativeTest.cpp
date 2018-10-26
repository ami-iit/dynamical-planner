/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Constraints.h>
#include <DynamicalPlannerPrivate/DynamicalConstraints.h>
#include <DynamicalPlannerPrivate/HyperbolicSecant.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/OptimalControlProblem.h>
#include <URDFdir.h>
#include <memory>
#include <string>
#include <cassert>

using namespace DynamicalPlanner::Private;

typedef struct {
    std::shared_ptr<DynamicalConstraints> dynamical;
    std::shared_ptr<CentroidalMomentumConstraint> centroidalMomentum;
    std::shared_ptr<CoMPositionConstraint> comPosition;
    std::vector<std::shared_ptr<NormalVelocityControlConstraints>> leftNormalVelocityControl, rightNormalVelocityControl;
    std::vector<std::shared_ptr<PlanarVelocityControlConstraints>> leftPlanarVelocityControl, rightPlanarVelocityControl;
    std::vector<std::shared_ptr<ContactForceControlConstraints>> leftContactsForceControl, rightContactsForceControl;
    std::vector<std::shared_ptr<ContactFrictionConstraint>> leftContactsFriction, rightContactsFriction;
    std::vector<std::shared_ptr<ContactPositionConsistencyConstraint>> leftContactsPosition, rightContactsPosition;
    std::shared_ptr<FeetLateralDistanceConstraint> feetLateralDistance;
    std::shared_ptr<QuaternionNormConstraint> quaternionNorm;
    std::vector<std::shared_ptr<DynamicalComplementarityConstraint>> leftComplementarity, rightComplementarity;
} ConstraintSet;

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

void setVariables(VariablesLabeller& stateVariables, VariablesLabeller& controlVariables, size_t numberOfPoints, size_t numberOfDofs) {

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

    ok = stateVariables.addLabel("JointsPosition", numberOfDofs);
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("BaseVelocity", 6);
    ASSERT_IS_TRUE(ok);

    ok = controlVariables.addLabel("JointsVelocity", numberOfDofs);
    ASSERT_IS_TRUE(ok);
}

void configureSharedKinDyn(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn) {
    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

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

void initializeConstraints(ConstraintSet& constraints, const std::vector<iDynTree::Position>& leftPositions,
                           const std::vector<iDynTree::Position>& rightPositions, const VariablesLabeller& stateVariables,
                           const VariablesLabeller& controlVariables, std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                           iDynTree::optimalcontrol::OptimalControlProblem& ocProblem) {
    iDynTree::Vector3 velocityMaximumDerivative;

    double forceMaximumDerivative = 10.0;
    double forceDissipationRatios = 10.0;
    double complementarityDissipation = 10.0;
    iDynTree::toEigen(velocityMaximumDerivative).setConstant(10.0);

    HyperbolicSecant forceActivation, velocityActivationZ;
    HyperbolicTangent velocityActivationXY;
    forceActivation.setScaling(1.0);
    velocityActivationXY.setScaling(1.0);
    velocityActivationZ.setScaling(1.0);

    iDynTree::FrameIndex leftFrame = timelySharedKinDyn->model().getFrameIndex("l_sole"),
            rightFrame = timelySharedKinDyn->model().getFrameIndex("r_sole");

    ASSERT_IS_TRUE(leftPositions.size() == rightPositions.size());
    size_t numberOfPoints = leftPositions.size();

    constraints.leftNormalVelocityControl.resize(numberOfPoints);
    constraints.leftPlanarVelocityControl.resize(numberOfPoints);
    constraints.leftContactsForceControl.resize(numberOfPoints);
    constraints.leftContactsFriction.resize(numberOfPoints);
    constraints.leftContactsPosition.resize(numberOfPoints);
    constraints.leftComplementarity.resize(numberOfPoints);

    constraints.rightNormalVelocityControl.resize(numberOfPoints);
    constraints.rightPlanarVelocityControl.resize(numberOfPoints);
    constraints.rightContactsForceControl.resize(numberOfPoints);
    constraints.rightContactsFriction.resize(numberOfPoints);
    constraints.rightContactsPosition.resize(numberOfPoints);
    constraints.rightComplementarity.resize(numberOfPoints);


    bool ok = false;
    constraints.dynamical = std::make_shared<DynamicalConstraints>(stateVariables, controlVariables, timelySharedKinDyn, velocityActivationXY);
    ok = ocProblem.setDynamicalSystemConstraint(constraints.dynamical);
    ASSERT_IS_TRUE(ok);

    for (size_t i = 0; i < numberOfPoints; ++i) {
        constraints.leftNormalVelocityControl[i] = std::make_shared<NormalVelocityControlConstraints>(stateVariables, controlVariables,
                                                                                                      "Left", i, velocityActivationZ,
                                                                                                      velocityMaximumDerivative(2));
        ok = ocProblem.addConstraint(constraints.leftNormalVelocityControl[i]);
        ASSERT_IS_TRUE(ok);


        constraints.leftPlanarVelocityControl[i] = std::make_shared<PlanarVelocityControlConstraints>(stateVariables, controlVariables,
                                                                                                      "Left", i,velocityActivationXY,
                                                                                                      velocityMaximumDerivative(0),
                                                                                                      velocityMaximumDerivative(1));
        ok = ocProblem.addConstraint(constraints.leftPlanarVelocityControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateVariables, controlVariables, "Left",
                                                                                                   i, forceActivation, forceMaximumDerivative,
                                                                                                   forceDissipationRatios,
                                                                                                   0.5);
        ok = ocProblem.addConstraint(constraints.leftContactsForceControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftContactsFriction[i] = std::make_shared<ContactFrictionConstraint>(stateVariables, controlVariables, "Left", i);
        ok = ocProblem.addConstraint(constraints.leftContactsFriction[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftContactsPosition[i] = std::make_shared<ContactPositionConsistencyConstraint>(stateVariables, controlVariables,
                                                                                                     timelySharedKinDyn, leftFrame, "Left",
                                                                                                     leftPositions[i], i);
        ok = ocProblem.addConstraint(constraints.leftContactsPosition[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftComplementarity[i] = std::make_shared<DynamicalComplementarityConstraint>(stateVariables, controlVariables,
                                                                                                  "Left", i, complementarityDissipation);
        ok = ocProblem.addConstraint(constraints.leftComplementarity[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightNormalVelocityControl[i] = std::make_shared<NormalVelocityControlConstraints>(stateVariables, controlVariables,
                                                                                                      "Right", i, velocityActivationZ,
                                                                                                      velocityMaximumDerivative(2));
        ok = ocProblem.addConstraint(constraints.rightNormalVelocityControl[i]);
        ASSERT_IS_TRUE(ok);


        constraints.rightPlanarVelocityControl[i] = std::make_shared<PlanarVelocityControlConstraints>(stateVariables, controlVariables,
                                                                                                      "Right", i,velocityActivationXY,
                                                                                                      velocityMaximumDerivative(0),
                                                                                                      velocityMaximumDerivative(1));
        ok = ocProblem.addConstraint(constraints.rightPlanarVelocityControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateVariables, controlVariables, "Right",
                                                                                                   i, forceActivation, forceMaximumDerivative,
                                                                                                   forceDissipationRatios,
                                                                                                   0.5);
        ok = ocProblem.addConstraint(constraints.rightContactsForceControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsFriction[i] = std::make_shared<ContactFrictionConstraint>(stateVariables, controlVariables, "Right", i);
        ok = ocProblem.addConstraint(constraints.rightContactsFriction[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsPosition[i] = std::make_shared<ContactPositionConsistencyConstraint>(stateVariables, controlVariables,
                                                                                                     timelySharedKinDyn, rightFrame, "Right",
                                                                                                     rightPositions[i], i);
        ok = ocProblem.addConstraint(constraints.rightContactsPosition[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightComplementarity[i] = std::make_shared<DynamicalComplementarityConstraint>(stateVariables, controlVariables,
                                                                                                  "Right", i, complementarityDissipation);
        ok = ocProblem.addConstraint(constraints.rightComplementarity[i]);
        ASSERT_IS_TRUE(ok);
    }
    constraints.centroidalMomentum = std::make_shared<CentroidalMomentumConstraint>(stateVariables, controlVariables, timelySharedKinDyn);
    ok = ocProblem.addConstraint(constraints.centroidalMomentum);
    ASSERT_IS_TRUE(ok);

    constraints.comPosition = std::make_shared<CoMPositionConstraint>(stateVariables, controlVariables, timelySharedKinDyn);
    ok = ocProblem.addConstraint(constraints.comPosition);
    ASSERT_IS_TRUE(ok);

    constraints.feetLateralDistance = std::make_shared<FeetLateralDistanceConstraint>(stateVariables, controlVariables, timelySharedKinDyn,
                                                                                      1, rightFrame, leftFrame);
    ok = ocProblem.addConstraint(constraints.feetLateralDistance);
    ASSERT_IS_TRUE(ok);

    constraints.quaternionNorm = std::make_shared<QuaternionNormConstraint>(stateVariables, controlVariables);
    ok = ocProblem.addConstraint(constraints.quaternionNorm);
    ASSERT_IS_TRUE(ok);
}

void checkDynamicalConstraintDerivative(double time, const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
                                        double perturbation, ConstraintSet& constraints) {
    iDynTree::VectorDynSize originalDynamics, perturbedDynamics, perturbedState, perturbedControl, firstOrderTaylor;
    iDynTree::MatrixDynSize stateJacobian, controlJacobian;
    originalDynamics.resize(originalStateVector.size());
    perturbedDynamics = originalDynamics;
    firstOrderTaylor = originalDynamics;
    bool ok = constraints.dynamical->setControlInput(originalControlVector);
    ASSERT_IS_TRUE(ok);
    ok = constraints.dynamical->dynamics(originalStateVector, time, originalDynamics);
    ASSERT_IS_TRUE(ok);
//    std::cerr << "Original dynamics:" << std::endl << originalDynamics.toString() << std::endl;


    ok = constraints.dynamical->dynamicsStateFirstDerivative(originalStateVector, time, stateJacobian);
    ASSERT_IS_TRUE(ok);

    ok = constraints.dynamical->dynamicsControlFirstDerivative(originalStateVector, time, controlJacobian);
    ASSERT_IS_TRUE(ok);

    for (unsigned int i = 0; i < originalStateVector.size(); ++i) {
        perturbedState = originalStateVector;
        perturbedState(i) = perturbedState(i) + perturbation;

        ok = constraints.dynamical->dynamics(perturbedState, 0.0, perturbedDynamics);
        ASSERT_IS_TRUE(ok);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalDynamics) +
                iDynTree::toEigen(stateJacobian) * (iDynTree::toEigen(perturbedState) - iDynTree::toEigen(originalStateVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedDynamics, firstOrderTaylor, perturbation/50);
    }

    for (unsigned int i = 0; i < originalControlVector.size(); ++i) {
        perturbedControl = originalControlVector;
        perturbedControl(i) = perturbedControl(i) + perturbation;

        ok = constraints.dynamical->setControlInput(perturbedControl);
        ASSERT_IS_TRUE(ok);
        ok = constraints.dynamical->dynamics(originalStateVector, 0.0, perturbedDynamics);
        ASSERT_IS_TRUE(ok);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalDynamics) +
                iDynTree::toEigen(controlJacobian) * (iDynTree::toEigen(perturbedControl) - iDynTree::toEigen(originalControlVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedDynamics, firstOrderTaylor, perturbation/50);
    }
}

void checkConstraintsDerivative(double time, const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
                                double perturbation, iDynTree::optimalcontrol::OptimalControlProblem &ocProblem) {

    iDynTree::VectorDynSize originalConstraints, perturbedConstraints, perturbedState, perturbedControl, firstOrderTaylor;
    iDynTree::MatrixDynSize stateJacobian, controlJacobian;


    originalConstraints.resize(ocProblem.getConstraintsDimension());
    perturbedConstraints = originalConstraints;
    firstOrderTaylor = originalConstraints;
    bool ok = false;

    ok = ocProblem.constraintsEvaluation(time, originalStateVector, originalControlVector, originalConstraints);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.constraintsJacobianWRTState(time, originalStateVector, originalControlVector, stateJacobian);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.constraintsJacobianWRTControl(time, originalStateVector, originalControlVector, controlJacobian);
    ASSERT_IS_TRUE(ok);

    for (unsigned int i = 0; i < originalStateVector.size(); ++i) {

//        std::cerr << "State: " << i << std::endl;
        perturbedState = originalStateVector;
        perturbedState(i) = perturbedState(i) + perturbation;
        ok = ocProblem.constraintsEvaluation(time, perturbedState, originalControlVector, perturbedConstraints);
        ASSERT_IS_TRUE(ok);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalConstraints) +
                iDynTree::toEigen(stateJacobian) * (iDynTree::toEigen(perturbedState) - iDynTree::toEigen(originalStateVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedConstraints, firstOrderTaylor, perturbation/10);
    }

    for (unsigned int i = 0; i < originalControlVector.size(); ++i) {
        perturbedControl = originalControlVector;
        perturbedControl(i) = perturbedControl(i) + perturbation;

        ok = ocProblem.constraintsEvaluation(time, originalStateVector, perturbedControl, perturbedConstraints);
        ASSERT_IS_TRUE(ok);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalConstraints) +
                iDynTree::toEigen(controlJacobian) * (iDynTree::toEigen(perturbedControl) - iDynTree::toEigen(originalControlVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedConstraints, firstOrderTaylor, perturbation/10);
    }
}


int main() {

    VariablesLabeller stateVariables, controlVariables;
    std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn = std::make_shared<TimelySharedKinDynComputations>();
    ConstraintSet constraints;
    iDynTree::optimalcontrol::OptimalControlProblem ocProblem;

    configureSharedKinDyn(timelySharedKinDyn);
    setVariables(stateVariables, controlVariables, 4, timelySharedKinDyn->model().getNrOfDOFs());
    std::vector<iDynTree::Position> leftPositions, rightPositions;
    leftPositions.resize(4);
    iDynTree::toEigen(leftPositions[0]) << 0.125, -0.04, 0.0;
    iDynTree::toEigen(leftPositions[1]) << 0.125,  0.04, 0.0;
    iDynTree::toEigen(leftPositions[2]) << -0.063,  0.04, 0.0;
    iDynTree::toEigen(leftPositions[3]) << 0.063, -0.04, 0.0;
    rightPositions.resize(4);
    iDynTree::toEigen(rightPositions[0]) << 0.125,  0.04, 0.0;
    iDynTree::toEigen(rightPositions[1]) << 0.125, -0.04, 0.0;
    iDynTree::toEigen(rightPositions[2]) << -0.063, -0.04, 0.0;
    iDynTree::toEigen(rightPositions[3]) << 0.063,  0.04, 0.0;

    initializeConstraints(constraints, leftPositions, rightPositions, stateVariables, controlVariables, timelySharedKinDyn, ocProblem);

    iDynTree::VectorDynSize stateVector, controlVector;
    stateVector.resize(static_cast<unsigned int>(stateVariables.size()));
    iDynTree::getRandomVector(stateVector);
//    std::cerr << "Original state:" << std::endl << stateVector.toString() << std::endl;
    controlVector.resize(static_cast<unsigned int>(controlVariables.size()));
    iDynTree::getRandomVector(controlVector);
//    std::cerr << "Original control:" << std::endl << controlVector.toString() << std::endl;


    checkDynamicalConstraintDerivative(0.0, stateVector, controlVector, 0.01, constraints);

    checkConstraintsDerivative(0.0, stateVector, controlVector, 0.001, ocProblem);

    checkDynamicalConstraintDerivative(1.0, stateVector, controlVector, 0.01, constraints);

    checkConstraintsDerivative(1.0, stateVector, controlVector, 0.001, ocProblem);

    return EXIT_SUCCESS;
}
