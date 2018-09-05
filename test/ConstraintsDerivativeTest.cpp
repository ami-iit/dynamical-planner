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
    std::vector<std::shared_ptr<ContactVelocityControlConstraints>> leftContactsVelocityControl, rightContactsVelocityControl;
    std::vector<std::shared_ptr<ContactForceControlConstraints>> leftContactsForceControl, rightContactsForceControl;
    std::vector<std::shared_ptr<ContactFrictionConstraint>> leftContactsFriction, rightContactsFriction;
    std::vector<std::shared_ptr<ContactPositionConsistencyConstraint>> leftContactsPosition, rightContactsPosition;
    std::shared_ptr<FeetLateralDistanceConstraint> feetLateralDistance;
    std::shared_ptr<QuaternionNormConstraint> quaternionNorm;
} ConstraintSet;

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

void initializeConstraints(ConstraintSet& constraints, const std::vector<iDynTree::Position>& leftPositions,
                           const std::vector<iDynTree::Position>& rightPositions, const VariablesLabeller& stateVariables,
                           const VariablesLabeller& controlVariables, std::shared_ptr<SharedKinDynComputation> sharedKinDyn,
                           iDynTree::optimalcontrol::OptimalControlProblem& ocProblem) {
    iDynTree::Vector3 forceMaximumDerivative, forceDissipationRatios, velocityMaximumDerivative, velocityDissipatioRatio;

    iDynTree::toEigen(forceMaximumDerivative).setConstant(10.0);
    iDynTree::toEigen(forceDissipationRatios).setConstant(10.0);
    iDynTree::toEigen(velocityMaximumDerivative).setConstant(10.0);
    iDynTree::toEigen(velocityDissipatioRatio).setConstant(10.0);

    HyperbolicSecant forceActivation, velocityActivationXY, velocityActivationZ;
    forceActivation.setScaling(1.0);
    velocityActivationXY.setScaling(1.0);
    velocityActivationZ.setScaling(1.0);

    iDynTree::FrameIndex leftFrame = sharedKinDyn->model().getFrameIndex("l_sole"),
            rightFrame = sharedKinDyn->model().getFrameIndex("r_sole");

    ASSERT_IS_TRUE(leftPositions.size() == rightPositions.size());
    size_t numberOfPoints = leftPositions.size();

    constraints.leftContactsVelocityControl.resize(numberOfPoints);
    constraints.leftContactsForceControl.resize(numberOfPoints);
    constraints.leftContactsFriction.resize(numberOfPoints);
    constraints.leftContactsPosition.resize(numberOfPoints);

    constraints.rightContactsVelocityControl.resize(numberOfPoints);
    constraints.rightContactsForceControl.resize(numberOfPoints);
    constraints.rightContactsFriction.resize(numberOfPoints);
    constraints.rightContactsPosition.resize(numberOfPoints);


    bool ok = false;
    constraints.dynamical = std::make_shared<DynamicalConstraints>(stateVariables, controlVariables, sharedKinDyn);
    ok = ocProblem.setDynamicalSystemConstraint(constraints.dynamical);
    ASSERT_IS_TRUE(ok);

    for (size_t i = 0; i < numberOfPoints; ++i) {
        constraints.leftContactsVelocityControl[i] = std::make_shared<ContactVelocityControlConstraints>(stateVariables, controlVariables,
                                                                                                         "Left", i, velocityActivationXY,
                                                                                                         velocityActivationZ,
                                                                                                         velocityMaximumDerivative,
                                                                                                         velocityDissipatioRatio);
        ok = ocProblem.addConstraint(constraints.leftContactsVelocityControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateVariables, controlVariables, "Left",
                                                                                                   i, forceActivation, forceMaximumDerivative,
                                                                                                   forceDissipationRatios);
        ok = ocProblem.addConstraint(constraints.leftContactsForceControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftContactsFriction[i] = std::make_shared<ContactFrictionConstraint>(stateVariables, controlVariables, "Left", i);
        ok = ocProblem.addConstraint(constraints.leftContactsFriction[i]);
        ASSERT_IS_TRUE(ok);

        constraints.leftContactsPosition[i] = std::make_shared<ContactPositionConsistencyConstraint>(stateVariables, controlVariables,
                                                                                                     sharedKinDyn, leftFrame, "Left",
                                                                                                     leftPositions[i], i);
        ok = ocProblem.addConstraint(constraints.leftContactsPosition[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsVelocityControl[i] = std::make_shared<ContactVelocityControlConstraints>(stateVariables, controlVariables,
                                                                                                         "Right", i, velocityActivationXY,
                                                                                                         velocityActivationZ,
                                                                                                         velocityMaximumDerivative,
                                                                                                         velocityDissipatioRatio);
        ok = ocProblem.addConstraint(constraints.rightContactsVelocityControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateVariables, controlVariables, "Right",
                                                                                                   i, forceActivation, forceMaximumDerivative,
                                                                                                   forceDissipationRatios);
        ok = ocProblem.addConstraint(constraints.rightContactsForceControl[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsFriction[i] = std::make_shared<ContactFrictionConstraint>(stateVariables, controlVariables, "Right", i);
        ok = ocProblem.addConstraint(constraints.rightContactsFriction[i]);
        ASSERT_IS_TRUE(ok);

        constraints.rightContactsPosition[i] = std::make_shared<ContactPositionConsistencyConstraint>(stateVariables, controlVariables,
                                                                                                     sharedKinDyn, rightFrame, "Right",
                                                                                                     rightPositions[i], i);
        ok = ocProblem.addConstraint(constraints.rightContactsPosition[i]);
        ASSERT_IS_TRUE(ok);
    }
    constraints.centroidalMomentum = std::make_shared<CentroidalMomentumConstraint>(stateVariables, controlVariables, sharedKinDyn);
    ok = ocProblem.addConstraint(constraints.centroidalMomentum);
    ASSERT_IS_TRUE(ok);

    constraints.comPosition = std::make_shared<CoMPositionConstraint>(stateVariables, controlVariables, sharedKinDyn);
    ok = ocProblem.addConstraint(constraints.comPosition);
    ASSERT_IS_TRUE(ok);

    constraints.feetLateralDistance = std::make_shared<FeetLateralDistanceConstraint>(stateVariables, controlVariables, sharedKinDyn,
                                                                                      1, rightFrame, leftFrame);
    ok = ocProblem.addConstraint(constraints.feetLateralDistance);
    ASSERT_IS_TRUE(ok);

    constraints.quaternionNorm = std::make_shared<QuaternionNormConstraint>(stateVariables, controlVariables);
    ok = ocProblem.addConstraint(constraints.quaternionNorm);
    ASSERT_IS_TRUE(ok);
}

void checkDynamicalConstraintDerivative(const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
                                        double perturbation, ConstraintSet& constraints) {
    iDynTree::VectorDynSize originalDynamics, perturbedDynamics, perturbedState, perturbedControl, firstOrderTaylor;
    iDynTree::MatrixDynSize stateJacobian, controlJacobian;
    originalDynamics.resize(originalStateVector.size());
    perturbedDynamics = originalDynamics;
    firstOrderTaylor = originalDynamics;
    bool ok = constraints.dynamical->setControlInput(originalControlVector);
    ASSERT_IS_TRUE(ok);
    ok = constraints.dynamical->dynamics(originalStateVector, 0.0, originalDynamics);
    ASSERT_IS_TRUE(ok);
//    std::cerr << "Original dynamics:" << std::endl << originalDynamics.toString() << std::endl;


    ok = constraints.dynamical->dynamicsStateFirstDerivative(originalStateVector, 0.0, stateJacobian);
    ASSERT_IS_TRUE(ok);

    ok = constraints.dynamical->dynamicsControlFirstDerivative(originalStateVector, 0.0, controlJacobian);
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

void checkConstraintsDerivative(const iDynTree::VectorDynSize& originalStateVector, const iDynTree::VectorDynSize& originalControlVector,
                                double perturbation, iDynTree::optimalcontrol::OptimalControlProblem &ocProblem) {

    iDynTree::VectorDynSize originalConstraints, perturbedConstraints, perturbedState, perturbedControl, firstOrderTaylor;
    iDynTree::MatrixDynSize stateJacobian, controlJacobian;


    originalConstraints.resize(ocProblem.getConstraintsDimension());
    perturbedConstraints = originalConstraints;
    firstOrderTaylor = originalConstraints;
    bool ok = false;

    ok = ocProblem.constraintsEvaluation(0.0, originalStateVector, originalControlVector, originalConstraints);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.constraintsJacobianWRTState(0.0, originalStateVector, originalControlVector, stateJacobian);
    ASSERT_IS_TRUE(ok);

    ok = ocProblem.constraintsJacobianWRTControl(0.0, originalStateVector, originalControlVector, controlJacobian);
    ASSERT_IS_TRUE(ok);

    for (unsigned int i = 0; i < originalStateVector.size(); ++i) {
        perturbedState = originalStateVector;
        perturbedState(i) = perturbedState(i) + perturbation;

        ok = ocProblem.constraintsEvaluation(0.0, perturbedState, originalControlVector, perturbedConstraints);
        ASSERT_IS_TRUE(ok);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalConstraints) +
                iDynTree::toEigen(stateJacobian) * (iDynTree::toEigen(perturbedState) - iDynTree::toEigen(originalStateVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedConstraints, firstOrderTaylor, perturbation/10);
    }

    for (unsigned int i = 0; i < originalControlVector.size(); ++i) {
        perturbedControl = originalControlVector;
        perturbedControl(i) = perturbedControl(i) + perturbation;

        ok = ocProblem.constraintsEvaluation(0.0, originalStateVector, perturbedControl, perturbedConstraints);
        ASSERT_IS_TRUE(ok);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalConstraints) +
                iDynTree::toEigen(controlJacobian) * (iDynTree::toEigen(perturbedControl) - iDynTree::toEigen(originalControlVector));
        ASSERT_EQUAL_VECTOR_TOL(perturbedConstraints, firstOrderTaylor, perturbation/10);
    }
}


int main() {

    VariablesLabeller stateVariables, controlVariables;
    std::shared_ptr<SharedKinDynComputation> sharedKinDyn = std::make_shared<SharedKinDynComputation>();
    ConstraintSet constraints;
    iDynTree::optimalcontrol::OptimalControlProblem ocProblem;

    setVariables(stateVariables, controlVariables, 4);
    configureSharedKinDyn(sharedKinDyn);
    std::vector<iDynTree::Position> leftPositions, rightPositions;
    leftPositions.resize(4);
    rightPositions.resize(4);
    leftPositions[0] = iDynTree::Position(0.125, -0.04, 0.0);
    leftPositions[1] = iDynTree::Position(0.125,  0.04, 0.0);
    leftPositions[2] = iDynTree::Position(-0.063,  0.04, 0.0);
    leftPositions[3] = iDynTree::Position( 0.063, -0.04, 0.0);
    rightPositions[0] = iDynTree::Position(0.125,  0.04, 0.0);
    rightPositions[1] = iDynTree::Position(0.125, -0.04, 0.0);
    rightPositions[2] = iDynTree::Position(-0.063, -0.04, 0.0);
    rightPositions[3] = iDynTree::Position( 0.063,  0.04, 0.0);

    initializeConstraints(constraints, leftPositions, rightPositions, stateVariables, controlVariables, sharedKinDyn, ocProblem);

    iDynTree::VectorDynSize stateVector, controlVector;
    stateVector.resize(static_cast<unsigned int>(stateVariables.size()));
    iDynTree::getRandomVector(stateVector);
//    std::cerr << "Original state:" << std::endl << stateVector.toString() << std::endl;
    controlVector.resize(static_cast<unsigned int>(controlVariables.size()));
    iDynTree::getRandomVector(controlVector);
//    std::cerr << "Original control:" << std::endl << controlVector.toString() << std::endl;


    checkDynamicalConstraintDerivative(stateVector, controlVector, 0.01, constraints);

    checkConstraintsDerivative(stateVector, controlVector, 0.001, ocProblem);

    return EXIT_SUCCESS;
}
