/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/SharedKinDynComputations.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <URDFdir.h>
#include <memory>

using namespace DynamicalPlanner::Private;

void configureSharedKinDyn(std::shared_ptr<SharedKinDynComputations> sharedKinDyn, iDynTree::LinkNetExternalWrenches &linkExtForces) {
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
    assert(sharedKinDyn);
    const iDynTree::Model& model = modelLoader.model();
    ok = sharedKinDyn->loadRobotModel(model);
    ASSERT_IS_TRUE(ok);
//    ASSERT_IS_TRUE(sharedKinDyn->model().getNrOfDOFs() == 23);
    linkExtForces.resize(model);
}

void validateStaticForcesDerivative(RobotState& robotState, std::shared_ptr<SharedKinDynComputations> sharedKinDyn,
                                    const iDynTree::LinkNetExternalWrenches &linkExtForces) {
    double perturbationValue = 1e-2;
    RobotState perturbedState = robotState;
    iDynTree::FreeFloatingGeneralizedTorques generalizedTorques;
    iDynTree::VectorDynSize jointTorques, jointTorquesPerturbed, firstOrderTaylor;
    iDynTree::MatrixDynSize jointsDerivative;
    bool ok = sharedKinDyn->getStaticForces(robotState, linkExtForces, generalizedTorques);
    ASSERT_IS_TRUE(ok);
    jointTorques = generalizedTorques.jointTorques();

    ok = sharedKinDyn->getStaticForcesJointsDerivative(robotState, linkExtForces, jointsDerivative);

//    std::cerr << "Jacobian:" << std::endl << jointsDerivative.toString() << std::endl;

    for (unsigned int i = 0; i < robotState.s.size(); ++i) {
//        std::cerr << "State: " << i << std::endl;
        perturbedState = robotState;
        perturbedState.s(i) = robotState.s(i) + perturbationValue;
        //perturbedState.world_T_base = iDynTree::getRandomTransform();

        ok = sharedKinDyn->getStaticForces(perturbedState, linkExtForces, generalizedTorques);
        ASSERT_IS_TRUE(ok);
        jointTorquesPerturbed = generalizedTorques.jointTorques();

        firstOrderTaylor = jointTorques;
        iDynTree::toEigen(firstOrderTaylor) += iDynTree::toEigen(jointsDerivative) * (iDynTree::toEigen(perturbedState.s) -
                                                                                               iDynTree::toEigen(robotState.s));
        ASSERT_EQUAL_VECTOR_TOL(jointTorquesPerturbed, firstOrderTaylor, perturbationValue/10);
    }
}

int main() {

    std::shared_ptr<SharedKinDynComputations> sharedKinDyn = std::make_shared<SharedKinDynComputations>();
    iDynTree::LinkNetExternalWrenches linkExtForces;
    configureSharedKinDyn(sharedKinDyn, linkExtForces);

    RobotState robotState;

    robotState.base_velocity = iDynTree::getRandomTwist();
    robotState.s.resize(static_cast<unsigned int>(sharedKinDyn->model().getNrOfDOFs()));
    iDynTree::getRandomVector(robotState.s);
    robotState.s_dot.resize(static_cast<unsigned int>(sharedKinDyn->model().getNrOfDOFs()));
    iDynTree::getRandomVector(robotState.s_dot);
    iDynTree::Transform randomTransform = iDynTree::getRandomTransform();
    robotState.base_position = randomTransform.getPosition();
    robotState.base_quaternion = randomTransform.getRotation().asQuaternion();
    iDynTree::LinkIndex leftLink = sharedKinDyn->model().getFrameLink(sharedKinDyn->model().getFrameIndex("l_sole"));
    ASSERT_IS_TRUE(sharedKinDyn->model().isValidLinkIndex(leftLink));

    iDynTree::LinkIndex rightLink = sharedKinDyn->model().getFrameLink(sharedKinDyn->model().getFrameIndex("r_sole"));
    ASSERT_IS_TRUE(sharedKinDyn->model().isValidLinkIndex(rightLink));

    iDynTree::Wrench genericWrench;
    genericWrench.zero();
    iDynTree::getRandomVector(genericWrench, -10, 10);
    linkExtForces(leftLink) = genericWrench;
    iDynTree::getRandomVector(genericWrench, -10, 10);
    linkExtForces(rightLink) = genericWrench;

    validateStaticForcesDerivative(robotState, sharedKinDyn, linkExtForces);

    return EXIT_SUCCESS;
}
