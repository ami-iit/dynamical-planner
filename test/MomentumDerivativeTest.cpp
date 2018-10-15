/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <URDFdir.h>
#include <memory>

using namespace DynamicalPlanner::Private;

void configureSharedKinDyn(std::shared_ptr<SharedKinDynComputations> sharedKinDyn) {
    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch"});

    iDynTree::ModelLoader modelLoader;
    bool ok = modelLoader.loadModelFromFile(getAbsModelPath("iCubGenova04.urdf"));
    ASSERT_IS_TRUE(ok);
    ok = modelLoader.loadReducedModelFromFullModel(modelLoader.model(), vectorList);
    ASSERT_IS_TRUE(ok);
    assert(sharedKinDyn);
    ok = sharedKinDyn->loadRobotModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
//    ASSERT_IS_TRUE(sharedKinDyn->model().getNrOfDOFs() == 23);
}

void validateMomentumDerivative(RobotState& robotState, std::shared_ptr<SharedKinDynComputations> sharedKinDyn) {
    double perturbationValue = 1e-2;
    RobotState perturbedState = robotState;
    iDynTree::SpatialMomentum originalMomentum, perturbedMomentum;
    iDynTree::Vector6 firstOrderTaylor;
    iDynTree::MatrixDynSize momentumPartialDerivative;
    originalMomentum = sharedKinDyn->getLinearAngularMomentum(robotState, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
    bool ok = sharedKinDyn->getLinearAngularMomentumJointsDerivative(robotState, momentumPartialDerivative);

    ASSERT_IS_TRUE(ok);


    for (unsigned int i = 0; i < robotState.s.size(); ++i) {
        perturbedState = robotState;
        perturbedState.s(i) = robotState.s(i) + perturbationValue;
        perturbedState.world_T_base = iDynTree::getRandomTransform();

        perturbedMomentum = sharedKinDyn->getLinearAngularMomentum(perturbedState,
                                                                   iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalMomentum);
        iDynTree::toEigen(firstOrderTaylor) += iDynTree::toEigen(momentumPartialDerivative) * (iDynTree::toEigen(perturbedState.s) -
                                                                                               iDynTree::toEigen(robotState.s));
        ASSERT_EQUAL_VECTOR_TOL(perturbedMomentum, firstOrderTaylor, perturbationValue/10);
    }
}

void validateVelocityDerivative(RobotState& robotState, std::shared_ptr<SharedKinDynComputations> sharedKinDyn) {
    double perturbationValue = 1e-2;
    RobotState perturbedState = robotState;
    iDynTree::Twist originalVelocity, perturbedVelocity;
    iDynTree::Vector6 firstOrderTaylor;
    iDynTree::MatrixDynSize velocityPartialDerivative;

    for(size_t l = 0; l < sharedKinDyn->model().getNrOfLinks(); ++l) {
        originalVelocity = sharedKinDyn->getFrameVel(robotState, static_cast<iDynTree::FrameIndex>(l), iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        bool ok = sharedKinDyn->getFrameVelJointsDerivative(robotState, static_cast<iDynTree::FrameIndex>(l), velocityPartialDerivative);
        ASSERT_IS_TRUE(ok);

        for (unsigned int i = 0; i < robotState.s.size(); ++i) {
            perturbedState = robotState;
            perturbedState.s(i) = robotState.s(i) + perturbationValue;

            perturbedVelocity = sharedKinDyn->getFrameVel(perturbedState, static_cast<iDynTree::FrameIndex>(l), iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

            iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalVelocity);
            iDynTree::toEigen(firstOrderTaylor) += iDynTree::toEigen(velocityPartialDerivative) * (iDynTree::toEigen(perturbedState.s) -
                                                                                                   iDynTree::toEigen(robotState.s));
            ASSERT_EQUAL_VECTOR_TOL(perturbedVelocity, firstOrderTaylor, perturbationValue/10);
        }
    }
}

int main() {

    std::shared_ptr<SharedKinDynComputations> sharedKinDyn = std::make_shared<SharedKinDynComputations>();
    configureSharedKinDyn(sharedKinDyn);

    RobotState robotState;

    robotState.base_velocity = iDynTree::getRandomTwist();
    robotState.s.resize(static_cast<unsigned int>(sharedKinDyn->model().getNrOfDOFs()));
    iDynTree::getRandomVector(robotState.s);
    robotState.s_dot.resize(static_cast<unsigned int>(sharedKinDyn->model().getNrOfDOFs()));
    iDynTree::getRandomVector(robotState.s_dot);
    robotState.world_T_base = iDynTree::getRandomTransform();

    validateVelocityDerivative(robotState, sharedKinDyn);
    validateMomentumDerivative(robotState, sharedKinDyn);

    return EXIT_SUCCESS;
}
