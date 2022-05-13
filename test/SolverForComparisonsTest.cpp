/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Solver.h>
#include <DynamicalPlanner/Visualizer.h>
#include <DynamicalPlanner/RectangularFoot.h>
#include <DynamicalPlanner/PositionReferenceGenerator.h>
#include <DynamicalPlanner/Logger.h>
#include <DynamicalPlanner/Utilities.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/Optimizers/IpoptInterface.h>
#include <iDynTree/Optimizers/WorhpInterface.h>
#include <iDynTree/Integrators/ForwardEuler.h>
#include <URDFdir.h>
#include <FolderPath.h>
#include <chrono>
#include <ctime>
#include <sstream>
#include <filesystem>
#include <iDynTree/Core/Utils.h>

#include "cmdline.h"

double maximumComplementarity(const DynamicalPlanner::State &state) {

    double maximum = 0;
    double complementarity;

    for (size_t i = 1; i < state.leftContactPointsState.size(); ++i) {
        complementarity = state.leftContactPointsState[i].pointForce(2) * state.leftContactPointsState[i].pointPosition(2);
        if (complementarity > maximum) {
            maximum = complementarity;
        }
    }

    for (size_t i = 1; i < state.rightContactPointsState.size(); ++i) {
        complementarity = state.rightContactPointsState[i].pointForce(2) * state.rightContactPointsState[i].pointPosition(2);
        if (complementarity > maximum) {
            maximum = complementarity;
        }
    }

    return maximum;
}

class ExternalOptions
{
public:
    std::string solver = "mumps";

    DynamicalPlanner::ComplementarityType complementarity = DynamicalPlanner::ComplementarityType::Dynamical;

    double velocity = 0.05;

    double normalForceDissipationRatio = 250.0;

    double normalForceHyperbolicSecantScaling = 5000.0;

    double complementarityDissipation = 20.0;

    double dynamicComplementarityUpperBound = 0.05;

    double classicalComplementarityTolerance= 0.004;

    ExternalOptions(const cmdline::parser &cmd)
    {
        solver = cmd.get<std::string>("solver");
        std::string complementarity_string = cmd.get<std::string>("complementarity");

        if (complementarity_string == "dynamical")
        {
            complementarity = DynamicalPlanner::ComplementarityType::Dynamical;
        }
        else if (complementarity_string == "classical")
        {
            complementarity = DynamicalPlanner::ComplementarityType::Classical;
        }
        else if (complementarity_string == "hyperbolic")
        {
            complementarity = DynamicalPlanner::ComplementarityType::HyperbolicSecantInequality;
        }
        else
        {
            std::cout << "CONFIGURATION ERROR: THe only available options for complementarity are dynamical, classical and hyperbolic." <<std::endl;
            m_valid = false;
        }


        velocity = cmd.get<double>("velocity");
        normalForceDissipationRatio = cmd.get<double>("normalForceDissipationRatio");
        normalForceHyperbolicSecantScaling = cmd.get<double>("normalForceHyperbolicSecantScaling");
        complementarityDissipation = cmd.get<double>("complementarityDissipation");
        dynamicComplementarityUpperBound = cmd.get<double>("dynamicComplementarityUpperBound");
        classicalComplementarityTolerance = cmd.get<double>("classicalComplementarityTolerance");
        m_valid = true;
    }

    bool isValid() const{
        return m_valid;
    }

private:
    bool m_valid{false};
};

void addOptions(cmdline::parser &cmd)
{
    cmd.add<std::string>("solver", 's',
                         "Linear solver to use with ipopt",
                         false,
                         "mumps");

    cmd.add<std::string>("complementarity", 'c',
                         "Type of complementairity condition",
                         false,
                         "dynamical");

    cmd.add<double>("velocity", 'v',
                    "Linear walking velocity",
                    false,
                    0.05);

    cmd.add<double>("normalForceDissipationRatio", 0,
                    "Normal Force Dissipation Ratio (one of the two parameters of the hyperbolic complementarity method).",
                    false,
                    250.0);

    cmd.add<double>("normalForceHyperbolicSecantScaling", 0,
                    "Normal Force Hyperbolic Secant Scaling (one of the two parameters of the hyperbolic complementarity method).",
                    false,
                    5000.0);

    cmd.add<double>("complementarityDissipation", 0,
                    "The rate of dissipation for the complementarity (one of the two parameters of the dynamical complementarity method).",
                    false,
                    20.0);

    cmd.add<double>("dynamicComplementarityUpperBound", 0,
                    "The upper-bound for the dynamic complementairty (one of the two parameters of the dynamical complementarity method).",
                    false,
                    0.05);

    cmd.add<double>("classicalComplementarityTolerance", 0,
                    "The upper bound for the classical complementarity.",
                    false,
                    0.004);
}

void printParsedOptions(const cmdline::parser &cmd)
{
    std::cout << "--------------------------------" <<std::endl;
    std::cout << "Parsed options:" << std::endl;
    std::cout << "  - solver: " << cmd.get<std::string>("solver") << std::endl;
    std::cout << "  - complementarity: " << cmd.get<std::string>("complementarity") << std::endl;
    std::cout << "  - velocity: " << cmd.get<double>("velocity") << std::endl;
    std::cout << "  - normalForceDissipationRatio: " << cmd.get<double>("normalForceDissipationRatio") << std::endl;
    std::cout << "  - normalForceHyperbolicSecantScaling: " << cmd.get<double>("normalForceHyperbolicSecantScaling") << std::endl;
    std::cout << "  - complementarityDissipation: " << cmd.get<double>("complementarityDissipation") << std::endl;
    std::cout << "  - dynamicComplementarityUpperBound: " << cmd.get<double>("dynamicComplementarityUpperBound") << std::endl;
    std::cout << "  - classicalComplementarityTolerance: " << cmd.get<double>("classicalComplementarityTolerance") << std::endl;
    std::cout << "--------------------------------" <<std::endl;
}


int main(int argc, char** argv) {

    cmdline::parser cmd;
    addOptions(cmd);
    cmd.parse_check(argc, argv);
    printParsedOptions(cmd);

    ExternalOptions externalOptions(cmd);
    ASSERT_IS_TRUE(externalOptions.isValid());

    DynamicalPlanner::Solver solver;

    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
//                                         "r_shoulder_pitch", "r_shoulder_roll", "l_hip_pitch", "l_hip_roll", "l_hip_yaw",
//                                         "l_knee", "l_ankle_pitch", "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw",
//                                         "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
//                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

//    std::vector<std::string> vectorList({"r_hip_pitch"});

    iDynTree::ModelLoader modelLoader;
    bool ok = modelLoader.loadModelFromFile(getAbsModelPath("iCubGenova04.urdf"));
    ASSERT_IS_TRUE(ok);
    ok = modelLoader.loadReducedModelFromFullModel(modelLoader.model(), vectorList);
    ASSERT_IS_TRUE(ok);

    DynamicalPlanner::SettingsStruct settingsStruct = DynamicalPlanner::Settings::Defaults(modelLoader.model());

    DynamicalPlanner::RectangularFoot leftFoot, rightFoot;

    double d = 0.09;
    double l = 0.19;

    iDynTree::Position topLeftPositionOfLeft(0.1265,  0.041, -0.015);
    ok = leftFoot.setFoot(l, d, topLeftPositionOfLeft);
    ASSERT_IS_TRUE(ok);

    ok = leftFoot.getPoints(iDynTree::Transform::Identity(), settingsStruct.leftPointsPosition);
    ASSERT_IS_TRUE(ok);

    iDynTree::Position topLeftPositionOfRight(0.1265,  0.041, -0.015);
    ok = rightFoot.setFoot(l, d, topLeftPositionOfRight);
    ASSERT_IS_TRUE(ok);

    ok = rightFoot.getPoints(iDynTree::Transform::Identity(), settingsStruct.rightPointsPosition);
    ASSERT_IS_TRUE(ok);

    DynamicalPlanner::State initialState;

    initialState.resize(vectorList.size(), settingsStruct.leftPointsPosition.size());

    iDynTree::VectorDynSize desiredInitialJoints(static_cast<unsigned int>(modelLoader.model().getNrOfDOFs()));

    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, -7, 22, 11, 30, -7, 22, 11, 30, 5.082, 0.406, -0.131,
                                              -45.249, -26.454, -0.351, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351;

//    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, -7, 22, -7, 22, 5.082, 0.406, -0.131, -45.249,
//                                             -26.454, -0.351, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351;

//    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351, 5.082,
//                                               0.406, -0.131, -45.249, -26.454, -0.351;

    iDynTree::toEigen(desiredInitialJoints) *= iDynTree::deg2rad(1.0);

    ok = DynamicalPlanner::Utilities::FillDefaultInitialState(settingsStruct, desiredInitialJoints, leftFoot, rightFoot, initialState);
    ASSERT_IS_TRUE(ok);
    ok = DynamicalPlanner::Utilities::SetMinContactPointToZero(settingsStruct, initialState);
    ASSERT_IS_TRUE(ok);

//    leftFoot.getNormalRatiosFromCoP(0.005, -0.01, settingsStruct.desiredLeftRatios);
//    rightFoot.getNormalRatiosFromCoP(0.005, 0.01, settingsStruct.desiredRightRatios);


    settingsStruct.desiredJointsTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredInitialJoints);

    for (auto& joint : settingsStruct.jointsVelocityLimits) {
        joint.first = -2;
        joint.second = 2;
    }

    iDynTree::toEigen(settingsStruct.jointsRegularizationWeights).topRows<3>().setConstant(0.1);
    iDynTree::toEigen(settingsStruct.jointsRegularizationWeights).segment<8>(3).setConstant(10.0);
    iDynTree::toEigen(settingsStruct.jointsRegularizationWeights).bottomRows<12>().setConstant(1.0);

//    iDynTree::toEigen(settingsStruct.jointsVelocityCostWeights).topRows<3>().setConstant(100.0);

    double torsoVelocityLimit = 0.6;
    settingsStruct.jointsVelocityLimits[0].first = -torsoVelocityLimit;
    settingsStruct.jointsVelocityLimits[0].second = torsoVelocityLimit;
    settingsStruct.jointsVelocityLimits[1].first = -torsoVelocityLimit;
    settingsStruct.jointsVelocityLimits[1].second = torsoVelocityLimit;
    settingsStruct.jointsVelocityLimits[2].first = -torsoVelocityLimit;
    settingsStruct.jointsVelocityLimits[2].second = torsoVelocityLimit;

    double armsVelocityLimit = 0.5;
    for (size_t i = 3; i < 11; ++i) {
        settingsStruct.jointsVelocityLimits[i].first = -armsVelocityLimit;
        settingsStruct.jointsVelocityLimits[i].second = armsVelocityLimit;
    }

    settingsStruct.jointsLimits[4].first = iDynTree::deg2rad(+10.0); //l_shoulder_roll
    settingsStruct.jointsLimits[8].first = iDynTree::deg2rad(+10.0); // r_shoulder_roll
//    settingsStruct.jointsLimits[12].first = iDynTree::deg2rad(-5.0); //l_hip_roll
//    settingsStruct.jointsLimits[18].first = iDynTree::deg2rad(-5.0); // r_hip_roll

    settingsStruct.frameCostActive = true;
    settingsStruct.staticTorquesCostActive = false;
    settingsStruct.forceMeanCostActive = true;
    settingsStruct.comCostActive = false;
    settingsStruct.comVelocityCostActive = true;
    settingsStruct.forceDerivativeCostActive = true;
    settingsStruct.pointAccelerationCostActive = true;
    settingsStruct.jointsRegularizationCostActive = false;
    settingsStruct.jointsVelocityCostActive = false;
    settingsStruct.swingCostActive = true;
    settingsStruct.phantomForcesCostActive = false;
    settingsStruct.meanPointPositionCostActive = true;
    settingsStruct.leftFootYawCostActive = true;
    settingsStruct.rightFootYawCostActive = true;
    settingsStruct.feetDistanceCostActive = false;
    settingsStruct.jointsVelocityForPosturalCostActive = true;
    settingsStruct.complementarityCostActive = false;
    settingsStruct.basePositionCostActive = false;
    settingsStruct.frameAngularVelocityCostActive = false;
    settingsStruct.baseQuaternionCostActive = true;
    settingsStruct.forceRatioCostActive = true;
    settingsStruct.baseLinearVelocityCostActive = false;
    settingsStruct.baseQuaternionVelocityCostActive = true;
    settingsStruct.angularMomentumCostActive = false;


    settingsStruct.frameCostOverallWeight = 90.0;
    settingsStruct.frameCostTimeVaryingWeight = std::make_shared<DynamicalPlanner::Utilities::ExponentialWeight>(100.0, -10.0, 1.0, 0.0);
    settingsStruct.jointsVelocityCostOverallWeight = 1e-1;
    settingsStruct.staticTorquesCostOverallWeight = 1e-5;
    settingsStruct.jointsRegularizationCostOverallWeight = 1e-1;
    settingsStruct.forceMeanCostOverallWeight = 1e-1; //5e-3 for exploiting more partial contacts
    settingsStruct.forceDerivativesCostOverallWeight = 1e-4;
    settingsStruct.pointAccelerationCostOverallWeight = 5.0;
    settingsStruct.pointAccelerationWeights(0) = 1.0;
    settingsStruct.pointAccelerationWeights(1) = 1.0;
    settingsStruct.pointAccelerationWeights(2) = 5.0;

    settingsStruct.swingCostOverallWeight = 1000;
    settingsStruct.swingCostWeights(0) = 1.0;
    settingsStruct.swingCostWeights(1) = 1.0;
    settingsStruct.swingCostWeights(2) = 1.0;
    settingsStruct.phantomForcesCostOverallWeight = 1.0;
    settingsStruct.meanPointPositionCostOverallWeight = 100;
    settingsStruct.comCostOverallWeight = 100;
    settingsStruct.comWeights(0) = 1.0;
    settingsStruct.comWeights(1) = 1.0;
    settingsStruct.comWeights(2) = 1.0;
    settingsStruct.comVelocityCostOverallWeight = 1.0;
    settingsStruct.comVelocityWeights(0) = 10.0;
    settingsStruct.comVelocityWeights(1) = 0.1;
    settingsStruct.comVelocityWeights(2) = 1.0;
    settingsStruct.leftFootYawCostOverallWeight = 2000.0;
    settingsStruct.rightFootYawCostOverallWeight = 2000.0;
    settingsStruct.feetDistanceCostOverallWeight = 1.0;
    settingsStruct.jointsVelocityForPosturalCostOverallWeight = 1e-1;
    settingsStruct.complementarityCostOverallWeight = 1e-3;
    settingsStruct.frameAngularVelocityCostOverallWeight = 0.1;//1.0;
    settingsStruct.rotationalPIDgain = 0.0;//10.0;
    settingsStruct.baseQuaternionCostOverallWeight = 50.0;
    settingsStruct.forceRatioCostOverallWeight = 1e-1;
    settingsStruct.baseLinearVelocityCostOverallWeight = 1.0;
    settingsStruct.baseLinearVelocityCostWeights(0) = 10.0;
    settingsStruct.baseLinearVelocityCostWeights(1) = 0.1;
    settingsStruct.baseLinearVelocityCostWeights(2) = 1.0;
    settingsStruct.baseQuaternionVelocityCostOverallWeight = 1e-3;
    settingsStruct.angularMomentumCostOverallWeight = 1e-3;

//    settingsStruct.minimumDt = 0.01;
//    settingsStruct.controlPeriod = 0.1;
//    settingsStruct.maximumDt = 0.05;

    settingsStruct.minimumDt = 0.1;
    settingsStruct.controlPeriod = 0.1;
    settingsStruct.maximumDt = 1.0;
    settingsStruct.horizon = 2.0;
    settingsStruct.activeControlPercentage = 1.0;

    settingsStruct.comCostActiveRange.setTimeInterval(settingsStruct.horizon*0, settingsStruct.horizon);
    //    auto comReference = std::make_shared<CoMReference>(initialState.comPosition, 0.2, 0.0, 0.0);
    iDynTree::VectorDynSize comPointReference(3);
    iDynTree::toEigen(comPointReference) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.1, 0.01, 0.0));
    auto comReference = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comPointReference);
    settingsStruct.desiredCoMTrajectory  = comReference;

    iDynTree::VectorDynSize comVelocityReference(3);
    iDynTree::toEigen(comVelocityReference) = iDynTree::toEigen(iDynTree::Position(externalOptions.velocity, 0.0, 0.0));
    auto comVelocityTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comVelocityReference);
    settingsStruct.desiredCoMVelocityTrajectory  = comVelocityTrajectory;

    iDynTree::VectorDynSize baseVelocityReference(3);
    iDynTree::toEigen(baseVelocityReference) = iDynTree::toEigen(iDynTree::Position(-0.05, 0.0, 0.0));//The velocity is expressed in base coordinates. The iCub root link has the x axis pointing backward
    auto baseVelocityTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(baseVelocityReference);
    settingsStruct.desiredBaseLinearVelocityTrajectory = baseVelocityTrajectory;

    settingsStruct.meanPointPositionCostActiveRange.setTimeInterval(settingsStruct.horizon * 0, settingsStruct.horizon);

    DynamicalPlanner::Utilities::SimpleWalkingStateMachine stateMachine;
    iDynTree::Position initialReference;
    iDynTree::toEigen(initialReference) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.1, 0.0, 0.0));
    initialReference(2) = 0.0;
    iDynTree::Position stepIncrement(0.1, 0.0, 0.0);
    ok = stateMachine.initialize(initialReference, stepIncrement,
                                 settingsStruct.horizon, settingsStruct.horizon,
                                 settingsStruct.minimumDt, 30.0, 30.0, 1.0, 40.0);
    ASSERT_IS_TRUE(ok);
    stateMachine.setVerbose(true);
    settingsStruct.desiredMeanPointPosition = stateMachine.references()->timeVaryingReference();
    settingsStruct.meanPointPositionCostTimeVaryingWeight = stateMachine.references()->timeVaryingWeight();

    settingsStruct.constrainTargetCoMPosition = false;
    settingsStruct.targetCoMPositionTolerance = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(0.02);
    settingsStruct.constrainTargetCoMPositionRange.setTimeInterval(settingsStruct.horizon * 0.6, settingsStruct.horizon);

    settingsStruct.comPositionConstraintTolerance = 1e-5*0;
    settingsStruct.centroidalMomentumConstraintTolerance = 1e-5*0;
    settingsStruct.quaternionModulusConstraintTolerance = 1e-2*0;
    settingsStruct.pointPositionConstraintTolerance = 1e-4*0;

    iDynTree::toEigen(settingsStruct.forceMaximumDerivative).setConstant(100.0);
//    settingsStruct.forceMaximumDerivative(0) = 10.0;
//    settingsStruct.forceMaximumDerivative(1) = 10.0;

    //ContactFrictionConstraint
    settingsStruct.frictionCoefficient = 0.3;

    settingsStruct.minimumFeetDistance = 0.10;
    settingsStruct.feetMaximumRelativeHeight = 0.04;
    settingsStruct.desiredSwingHeight = 0.02;

    settingsStruct.maximumAngularMomentum = 10;


    //ContactVelocityControlConstraints
    iDynTree::toEigen(settingsStruct.velocityMaximumDerivative).setConstant(5.0);
    settingsStruct.velocityMaximumDerivative(0) = 2.0;
    settingsStruct.velocityMaximumDerivative(1) = 2.0;
    settingsStruct.planarVelocityHyperbolicTangentScaling = 10.0; //scales the position along z
    settingsStruct.normalVelocityHyperbolicSecantScaling = 5.0; //scales the force along z
    settingsStruct.classicalPlanarComplementarityTolerance = 0.02;
    settingsStruct.planarComplementarity = DynamicalPlanner::PlanarComplementarityType::HyperbolicTangentInDynamics;

    if (settingsStruct.planarComplementarity == DynamicalPlanner::PlanarComplementarityType::HyperbolicTangentInequality)
    {
        settingsStruct.velocityMaximumDerivative(0) = 4.0;
        settingsStruct.velocityMaximumDerivative(1) = 4.0;
        settingsStruct.planarVelocityHyperbolicTangentScaling = 5.0; //scales the position along z
    }

    settingsStruct.complementarity = externalOptions.complementarity;
    settingsStruct.normalForceDissipationRatio = externalOptions.normalForceDissipationRatio;
    settingsStruct.normalForceHyperbolicSecantScaling = externalOptions.normalForceHyperbolicSecantScaling;
    settingsStruct.complementarityDissipation = externalOptions.complementarityDissipation;
    settingsStruct.dynamicComplementarityUpperBound = externalOptions.dynamicComplementarityUpperBound;
    settingsStruct.classicalComplementarityTolerance = externalOptions.classicalComplementarityTolerance;

    settingsStruct.minimumCoMHeight = 0.5 * initialState.comPosition(2);

    iDynTree::VectorDynSize desiredQuaternion(4);
    desiredQuaternion = initialState.worldToBaseTransform.getRotation().asQuaternion();
    settingsStruct.desiredBaseQuaternionTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredQuaternion);

    auto ipoptSolver = std::make_shared<iDynTree::optimization::IpoptInterface>();

    ASSERT_IS_TRUE(ipoptSolver->isAvailable());

    ok = ipoptSolver->setIpoptOption("linear_solver", externalOptions.solver);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("ma57_pivtol", 1e-6);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("print_level", 5);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("nlp_scaling_max_gradient", 100.0);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("nlp_scaling_min_value", 1e-6);
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("evaluate_orig_obj_at_resto_trial", "no");
//    ASSERT_IS_TRUE(ok);

//    ok = ipoptSolver->setIpoptOption("check_derivatives_for_naninf", "yes");
//    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("tol", 1e-3);
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("bound_relax_factor", 1e-5);
//    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("honor_original_bounds", "yes");
//    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("dual_inf_tol", 1000.0);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("compl_inf_tol", 1e-2);
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("recalc_y","yes");
//    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("constr_viol_tol", 1e-4);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("acceptable_tol", 1e0);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("acceptable_iter", 2);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("acceptable_compl_inf_tol", 1.0);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("alpha_for_y", "dual-and-full");
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("accept_every_trial_step", "yes");
//    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("max_iter", 4000);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("ma97_print_level", -10);
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("mu_init", 1e3);
//    ASSERT_IS_TRUE(ok);

    ok = ipoptSolver->setIpoptOption("warm_start_bound_frac", 1e-2);
    ASSERT_IS_TRUE(ok);

    ok = ipoptSolver->setIpoptOption("warm_start_bound_push", 1e-2);
    ASSERT_IS_TRUE(ok);

    ok = ipoptSolver->setIpoptOption("warm_start_mult_bound_push", 1e-2);
    ok = ipoptSolver->setIpoptOption("warm_start_slack_bound_frac", 1e-2);
    ok = ipoptSolver->setIpoptOption("warm_start_slack_bound_push", 1e-2);
    ok = ipoptSolver->setIpoptOption("warm_start_init_point", "yes");
//    ok = ipoptSolver->setIpoptOption("warm_start_same_structure", "no");
//    ok = ipoptSolver->setIpoptOption("expect_infeasible_problem", "yes");
    ok = ipoptSolver->setIpoptOption("required_infeasibility_reduction", 0.8);

//    ok = ipoptSolver->setIpoptOption("first_hessian_perturbation", 100.0);
    ok = ipoptSolver->setIpoptOption("perturb_dec_fact", 0.1);
    ok = ipoptSolver->setIpoptOption("max_hessian_perturbation", 100.0);
//    ok = ipoptSolver->setIpoptOption("jacobian_regularization_value", 0.001);
//    ok = ipoptSolver->setIpoptOption("obj_scaling_factor", 500.0);

    ok = ipoptSolver->setIpoptOption("fast_step_computation", "yes");
//    ok = ipoptSolver->setIpoptOption("evaluate_orig_obj_at_resto_trial", "yes");
//    ok = ipoptSolver->setIpoptOption("soft_resto_pderror_reduction_factor", 0.0);
//    ok = ipoptSolver->setIpoptOption("linear_system_scaling", "slack-based");

    ipoptSolver->useApproximatedHessians(true);

    ok = solver.setOptimizer(ipoptSolver);
    ASSERT_IS_TRUE(ok);

    ok = solver.specifySettings(settingsStruct);
    ASSERT_IS_TRUE(ok);

    ok = solver.setInitialState(initialState);
    ASSERT_IS_TRUE(ok);

    DynamicalPlanner::State initialStateForGuess(initialState);

    for (auto& point : initialStateForGuess.leftContactPointsState)
    {
        iDynTree::toEigen(point.pointPosition) += 2 * iDynTree::toEigen(stepIncrement);
        point.pointForce.zero();

    }
    iDynTree::toEigen(initialStateForGuess.comPosition) += 0.1 * iDynTree::toEigen(stepIncrement);
//    initialStateForGuess.comPosition(1) -= std::abs(initialState.rightContactPointsState.begin()->pointPosition(1) - initialState.leftContactPointsState.begin()->pointPosition(1))/2.0;

    auto stateGuesses = std::make_shared<DynamicalPlanner::TimeInvariantState>(initialStateForGuess);

    DynamicalPlanner::Control initialControl(vectorList.size(), settingsStruct.leftPointsPosition.size());

    auto controlGuesses = std::make_shared<DynamicalPlanner::TimeInvariantControl>(initialControl);

    ok = solver.setGuesses(stateGuesses, controlGuesses);
    ASSERT_IS_TRUE(ok);

    std::vector<DynamicalPlanner::State> optimalStates;
    std::vector<DynamicalPlanner::Control> optimalControls;

    DynamicalPlanner::Visualizer visualizer, instanceVisualizer;
    ok = visualizer.setModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
    visualizer.setCameraPosition(iDynTree::Position(1.5, 0.0, 0.0));
    ok = visualizer.visualizeState(initialState);
    ASSERT_IS_TRUE(ok);

    ok = instanceVisualizer.setModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
    instanceVisualizer.setCameraPosition(iDynTree::Position(1.5, 0.0, 0.0));
    ok = instanceVisualizer.visualizeState(initialState);
    instanceVisualizer.visualizeWorldFrame(false);
    ASSERT_IS_TRUE(ok);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    ok = solver.solve(optimalStates, optimalControls);
    ASSERT_IS_TRUE(ok);
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time (1st): " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0 <<std::endl;

    auto timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm timeStruct;
    timeStruct = *std::localtime(&timeNow);
    std::ostringstream timeString;
    timeString << timeStruct.tm_year + 1900 << "-" << timeStruct.tm_mon + 1 << "-";
    timeString << timeStruct.tm_mday << "_" << timeStruct.tm_hour << "_" << timeStruct.tm_min;
    timeString << "_" << timeStruct.tm_sec;

    std::ostringstream complementarityString;
    switch (settingsStruct.complementarity)
    {
    case DynamicalPlanner::ComplementarityType::Classical:
        complementarityString.precision(3);
        complementarityString << "Classical_" << std::fixed << settingsStruct.classicalComplementarityTolerance;
        break;
    case DynamicalPlanner::ComplementarityType::Dynamical:
        complementarityString.precision(2);
        complementarityString << "Dynamical_ub" << std::fixed << settingsStruct.dynamicComplementarityUpperBound;
        complementarityString.precision(1);
        complementarityString << "_K" << std::fixed << settingsStruct.complementarityDissipation;
        break;
    case DynamicalPlanner::ComplementarityType::HyperbolicSecantInDynamics:
        complementarityString.precision(1);
        complementarityString << "HyperbolicSecantInDynamics_sec" << std::fixed << settingsStruct.normalForceHyperbolicSecantScaling
                              << "_K" << std::fixed << settingsStruct.normalForceDissipationRatio;
        break;
    case DynamicalPlanner::ComplementarityType::HyperbolicSecantInequality:
        complementarityString.precision(1);
        complementarityString << "HyperbolicSecantInequality_sec" << std::fixed << settingsStruct.normalForceHyperbolicSecantScaling
                              << "_K" << std::fixed << settingsStruct.normalForceDissipationRatio;
        break;
    }

    std::ostringstream velocityString;
    velocityString.precision(3);
    velocityString << "incr-"  << std::fixed << iDynTree::toEigen(stepIncrement).norm()
                   << "_speed-" << std::fixed << iDynTree::toEigen(comVelocityReference).norm();

    std::string outputFolder = getAbsDirPath("SavedVideos") + "/" + timeString.str() + "_" + complementarityString.str() + "_" + velocityString.str();
    std::filesystem::create_directories(outputFolder);
    ok = visualizer.visualizeStatesAndSaveAnimation(optimalStates, outputFolder, "test-1stIteration-" + timeString.str(), "mp4", settingsStruct.horizon * settingsStruct.activeControlPercentage);
    ASSERT_IS_TRUE(ok);

    //-----------------------

    ok = ipoptSolver->setIpoptOption("print_level", 0);
    ASSERT_IS_TRUE(ok);

    std::vector<DynamicalPlanner::State> mpcStates;
    std::vector<DynamicalPlanner::Control> mpcControls;
    std::vector<std::vector<DynamicalPlanner::State>> fullStates;
    std::vector<double> durations;
    std::vector<double> costs;
    mpcStates.push_back(initialState);

    visualizer.setCameraPosition(iDynTree::Position(2.0, 0.5, 0.5));
    double runningMean = 0;
    double runningCostMean = 0;
    double currentDuration;
    for (size_t i = 0; i < 100; ++i) {
        double initialTime;
        initialState = mpcStates.back();
        initialTime = initialState.time;
        initialState.time = 0.0;
        ok = DynamicalPlanner::Utilities::SetMinContactPointToZero(settingsStruct, initialState);
        ASSERT_IS_TRUE(ok);
        ok = solver.setInitialState(initialState);
        ASSERT_IS_TRUE(ok);
        begin = std::chrono::steady_clock::now();
        ok = solver.solve(optimalStates, optimalControls);
        if (!ok)
            break;
        end= std::chrono::steady_clock::now();
        currentDuration = (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0;
        runningMean = (runningMean * i + currentDuration)/(i+1);
        std::cout << "Elapsed time (" << i << "): " << currentDuration <<std::endl;
        std::cout << "Complementarity: " << maximumComplementarity(optimalStates.front()) << std::endl;
        std::cout << "Mean Time: " << runningMean << std::endl;
        double optimalCost(-1.0);
        ipoptSolver->getOptimalCost(optimalCost);
        runningCostMean = (runningCostMean * i + optimalCost) / (i+1);
        std::cout << "Cost Value: " << optimalCost << " (Mean: " << runningCostMean << ")" << std::endl;

        durations.push_back(currentDuration);
        costs.push_back(optimalCost);
        mpcStates.push_back(optimalStates.front());
        mpcStates.back().time += initialTime;
        mpcStates.push_back(optimalStates[1]);
        mpcStates.back().time += initialTime;
        mpcControls.push_back(optimalControls.front());
        mpcControls.back().time += initialTime;
        mpcControls.push_back(optimalControls[1]);
        mpcControls.back().time += initialTime;
        fullStates.push_back(optimalStates);
        iDynTree::Position instanceCameraPosition;
        iDynTree::toEigen(instanceCameraPosition) = iDynTree::toEigen(mpcStates.back().comPosition) + Eigen::Vector3d(1.5, 0.0, -mpcStates.back().comPosition(2));
        instanceVisualizer.setCameraPosition(instanceCameraPosition);
        instanceVisualizer.visualizeStates(optimalStates);
        visualizer.visualizeState(mpcStates.back());

        ok = stateMachine.advance(optimalStates.front(), optimalStates[static_cast<size_t>(std::round(optimalStates.size() * 0.3))]);
        ASSERT_IS_TRUE(ok);
    }

    timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    timeStruct = *std::localtime(&timeNow);
    timeString << timeStruct.tm_year + 1900 << "-" << timeStruct.tm_mon + 1 << "-";
    timeString << timeStruct.tm_mday << "_" << timeStruct.tm_hour << "_" << timeStruct.tm_min;
    timeString << "_" << timeStruct.tm_sec;

    ok = visualizer.visualizeStatesAndSaveAnimation(mpcStates, outputFolder, "test-" + timeString.str(), "mp4", settingsStruct.horizon * settingsStruct.activeControlPercentage);
    ASSERT_IS_TRUE(ok);

    auto stateCameraControl = [](const DynamicalPlanner::State&) {return iDynTree::Position(2.0, 0.5, 0.5);};
    auto instanceCameraControl = [](const DynamicalPlanner::State& inputState) {
        iDynTree::Position instanceCameraPosition;
        iDynTree::toEigen(instanceCameraPosition) = iDynTree::toEigen(inputState.comPosition) + Eigen::Vector3d(1.5, 0.0, -inputState.comPosition(2));
        return instanceCameraPosition;
    };

    visualizer.visualizeWorldFrame(false);
    ok = visualizer.visualizeMPCStatesAndSaveAnimation(mpcStates, stateCameraControl, fullStates, instanceCameraControl,
                                                       outputFolder, "test-SC-" + timeString.str(), "mp4", 2.0, settingsStruct.horizon * settingsStruct.activeControlPercentage);
    ASSERT_IS_TRUE(ok);

    DynamicalPlanner::Logger::saveSolutionVectorsToFile(outputFolder + "/log-" + timeString.str() + ".mat" , settingsStruct, mpcStates, mpcControls, durations, costs);

    return EXIT_SUCCESS;
}
