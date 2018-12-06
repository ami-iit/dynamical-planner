/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Solver.h>
#include <DynamicalPlanner/Visualizer.h>
#include <DynamicalPlanner/RectangularFoot.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/KinDynComputations.h>
#include <iDynTree/Optimizers/IpoptInterface.h>
#include <iDynTree/Optimizers/WorhpInterface.h>
#include <iDynTree/Integrators/ForwardEuler.h>
#include <URDFdir.h>
#include <FolderPath.h>
#include <chrono>
#include <ctime>
#include <sstream>

void fillInitialState(iDynTree::KinDynComputations& kinDyn, const DynamicalPlanner::SettingsStruct& settings,
                      const iDynTree::VectorDynSize& desiredJoints, DynamicalPlanner::RectangularFoot &foot,
                      DynamicalPlanner::State &initialState) {

    const iDynTree::Model& model = kinDyn.model();

    bool ok = kinDyn.setFloatingBase(model.getLinkName(model.getFrameLink(model.getFrameIndex(settings.leftFrameName))));
    ASSERT_IS_TRUE(ok);

    iDynTree::Vector3 gravity;
    gravity.zero();
    gravity(2) = -9.81;

    iDynTree::Twist baseVelocity = iDynTree::Twist::Zero();
    baseVelocity(0) = 0.3*0;

    ok = kinDyn.setRobotState(model.getFrameTransform(model.getFrameIndex(settings.leftFrameName)).inverse(), desiredJoints,
                              baseVelocity, iDynTree::VectorDynSize(desiredJoints.size()), gravity);
    ASSERT_IS_TRUE(ok);

    initialState.comPosition = kinDyn.getCenterOfMassPosition();

    initialState.jointsConfiguration = desiredJoints;

    kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
    initialState.momentumInCoM = kinDyn.getCentroidalTotalMomentum().asVector();

    initialState.time = 0.0;

    initialState.worldToBaseTransform = kinDyn.getWorldTransform(settings.floatingBaseName);

    iDynTree::Transform leftTransform = kinDyn.getWorldTransform(settings.leftFrameName);
    iDynTree::Transform rightTransform = kinDyn.getWorldTransform(settings.rightFrameName);

    double totalMass = 0.0;

    for(size_t l=0; l < model.getNrOfLinks(); l++)
    {
        totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
    }

    double normalForce = totalMass * 9.81;

    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointPosition = leftTransform * settings.leftPointsPosition[i];
        initialState.rightContactPointsState[i].pointPosition = rightTransform * settings.rightPointsPosition[i];
    }

    iDynTree::Wrench leftWrench, rightWrench;
    leftWrench.zero();
    rightWrench.zero();

    leftWrench(2) = normalForce/(1 + std::fabs(initialState.comPosition(1) - leftTransform.getPosition()(1))/std::fabs(initialState.comPosition(1) - rightTransform.getPosition()(1)));

    rightWrench(2) = normalForce - leftWrench(2);

    leftWrench(4) = -leftWrench(2) * (initialState.comPosition(0) - leftTransform.getPosition()(0));
    rightWrench(4) = -rightWrench(2) * (initialState.comPosition(0) - rightTransform.getPosition()(0));

    std::vector<iDynTree::Force> leftPointForces, rightPointForces;

    ok = foot.getForces(leftWrench, leftPointForces);
    ASSERT_IS_TRUE(ok);

    ok = foot.getForces(rightWrench, rightPointForces);
    ASSERT_IS_TRUE(ok);

    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointForce = leftPointForces[i];
        initialState.rightContactPointsState[i].pointForce =rightPointForces[i];
    }

}

void reconstructState(iDynTree::KinDynComputations& kinDyn, const DynamicalPlanner::SettingsStruct& settings, DynamicalPlanner::State &initialState) {

    bool ok = kinDyn.setFloatingBase(settings.floatingBaseName);
    ASSERT_IS_TRUE(ok);

    iDynTree::Vector3 gravity;
    gravity.zero();
    gravity(2) = -9.81;

    ok = kinDyn.setRobotState(initialState.worldToBaseTransform, initialState.jointsConfiguration,
                                   iDynTree::Twist::Zero(), iDynTree::VectorDynSize(initialState.jointsConfiguration.size()), gravity);
    ASSERT_IS_TRUE(ok);

    initialState.comPosition = kinDyn.getCenterOfMassPosition();

    initialState.time = 0.0;

    iDynTree::Transform leftTransform = kinDyn.getWorldTransform(settings.leftFrameName);
    iDynTree::Transform rightTransform = kinDyn.getWorldTransform(settings.rightFrameName);


    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointPosition = leftTransform * settings.leftPointsPosition[i];
        initialState.rightContactPointsState[i].pointPosition = rightTransform * settings.rightPointsPosition[i];
    }
}

class CoMReference : public iDynTree::optimalcontrol::TimeVaryingVector {
    iDynTree::VectorDynSize desiredCoM;
    double xVelocity, yVelocity, zVelocity;
    iDynTree::Vector3 initialCoM;

public:
    CoMReference(iDynTree::Vector3 &CoMinitial, double velX, double velY, double velZ)
        : desiredCoM(3)
        , xVelocity(velX)
        , yVelocity(velY)
        , zVelocity(velZ)
        , initialCoM(CoMinitial)
    { }

    ~CoMReference() override;

    iDynTree::VectorDynSize &get(double time, bool &isValid) override {
        desiredCoM(0) = initialCoM(0) + xVelocity * time;
        desiredCoM(1) = initialCoM(1) + yVelocity * time;
        desiredCoM(2) = initialCoM(2) + zVelocity * time;
        isValid = true;
        return desiredCoM;
    }

};
CoMReference::~CoMReference() { }

class StateGuess : public DynamicalPlanner::TimeVaryingState {
    DynamicalPlanner::State m_state, m_initialState;
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> m_comReference;
public:

    StateGuess(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> comReference, const DynamicalPlanner::State &initialState)
        : m_state(initialState)
        , m_initialState(initialState)
        , m_comReference(comReference)
    { }

    ~StateGuess() override;

    DynamicalPlanner::State &get(double time, bool &isValid) override {
        iDynTree::toEigen(m_state.comPosition) = iDynTree::toEigen(m_comReference->get(time, isValid));
        m_state.jointsConfiguration = m_initialState.jointsConfiguration;
        m_state.momentumInCoM.zero();
        m_state.worldToBaseTransform.setRotation(m_initialState.worldToBaseTransform.getRotation());
        iDynTree::Position basePosition, comDifference;
        iDynTree::toEigen(comDifference) = iDynTree::toEigen(m_state.comPosition) - iDynTree::toEigen(m_initialState.comPosition);
        iDynTree::toEigen(basePosition) = iDynTree::toEigen(m_initialState.worldToBaseTransform.getPosition()) + iDynTree::toEigen(comDifference);
        m_state.worldToBaseTransform.setPosition(basePosition);

        for (size_t i = 0; i < m_state.leftContactPointsState.size(); ++i) {
            iDynTree::toEigen(m_state.leftContactPointsState[i].pointPosition) = iDynTree::toEigen(m_initialState.leftContactPointsState[i].pointPosition) + iDynTree::toEigen(comDifference);
            m_state.leftContactPointsState[i].pointPosition(2) = 0;
        }

        for (size_t i = 0; i < m_state.rightContactPointsState.size(); ++i) {
            iDynTree::toEigen(m_state.rightContactPointsState[i].pointPosition) = iDynTree::toEigen(m_initialState.rightContactPointsState[i].pointPosition) + iDynTree::toEigen(comDifference);
            m_state.rightContactPointsState[i].pointPosition(2) = 0;
        }

        isValid = true;

        m_state.time = time;
        return m_state;
    }
};
StateGuess::~StateGuess() {}

int main() {

    DynamicalPlanner::Solver solver;
    DynamicalPlanner::Settings settings;

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

    DynamicalPlanner::RectangularFoot foot;

    double d = 0.08;
    double l = 0.188;

    iDynTree::Position topLeftPosition(0.125,  0.04, 0.0);
    ok = foot.setFoot(l, d, topLeftPosition);
    ASSERT_IS_TRUE(ok);

    ok = foot.getPoints(iDynTree::Transform::Identity(), settingsStruct.leftPointsPosition);
    ASSERT_IS_TRUE(ok);

    settingsStruct.rightPointsPosition = settingsStruct.leftPointsPosition;

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

    iDynTree::KinDynComputations kinDyn;

    ok = kinDyn.loadRobotModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);

    fillInitialState(kinDyn, settingsStruct, desiredInitialJoints, foot, initialState);
    reconstructState(kinDyn, settingsStruct, initialState);


    settingsStruct.desiredJointsTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredInitialJoints);

    for (auto& joint : settingsStruct.jointsVelocityLimits) {
        joint.first = -10;
        joint.second = 10;
    }

    settingsStruct.frameCostActive = true;
    settingsStruct.staticTorquesCostActive = true;
    settingsStruct.forceMeanCostActive = true;
    settingsStruct.comCostActive = false;
    settingsStruct.forceDerivativeCostActive = false;
    settingsStruct.pointAccelerationCostActive = false;
    settingsStruct.jointsRegularizationCostActive = true;
    settingsStruct.jointsVelocityCostActive = true;
    settingsStruct.swingCostActive = true;
    settingsStruct.phantomForcesCostActive = false;
    settingsStruct.meanPointPositionCostActive = true;

    settingsStruct.frameCostOverallWeight = 10;
    settingsStruct.jointsVelocityCostOverallWeight = 1e-2;
    settingsStruct.staticTorquesCostOverallWeight = 1e-5;
    settingsStruct.jointsRegularizationCostOverallWeight = 1e-1;
    settingsStruct.forceMeanCostOverallWeight = 1e-3;
    settingsStruct.forceDerivativesCostOverallWeight = 1e-15;
    settingsStruct.pointAccelerationCostOverallWeight = 1e-15;
    settingsStruct.swingCostOverallWeight = 10;
    settingsStruct.phantomForcesCostOverallWeight = 1.0;
    settingsStruct.meanPointPositionCostOverallWeight = 100;
    settingsStruct.comCostOverallWeight = 100;
    settingsStruct.comWeights(0) = 1.0;
    settingsStruct.comWeights(1) = 1.0;
    settingsStruct.comWeights(2) = 1.0;

//    settingsStruct.minimumDt = 0.01;
//    settingsStruct.controlPeriod = 0.1;
//    settingsStruct.maximumDt = 0.05;

    settingsStruct.minimumDt = 0.1;
    settingsStruct.controlPeriod = 0.1;
    settingsStruct.maximumDt = 1.0;
    settingsStruct.horizon = 1.5;
    settingsStruct.activeControlPercentage = 1.0;

    settingsStruct.comCostActiveRange.setTimeInterval(settingsStruct.horizon*0, settingsStruct.horizon);
    //    auto comReference = std::make_shared<CoMReference>(initialState.comPosition, 0.2, 0.0, 0.0);
    iDynTree::VectorDynSize comPointReference(3);
    iDynTree::toEigen(comPointReference) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.1, 0.01, 0.0));
    auto comReference = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comPointReference);
    settingsStruct.desiredCoMTrajectory  = comReference;

    settingsStruct.meanPointPositionCostActiveRange.setTimeInterval(settingsStruct.horizon * 0, settingsStruct.horizon);
    iDynTree::Position meanPointReference;
    iDynTree::toEigen(meanPointReference) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.2, 0.0, 0.0));
    meanPointReference(2) = 0.0;
    auto meanPointReferencePointer = std::make_shared<iDynTree::optimalcontrol::TimeInvariantPosition>(meanPointReference);
    settingsStruct.desiredMeanPointPosition = meanPointReferencePointer;

    settingsStruct.constrainTargetCoMPosition = false;
    settingsStruct.targetCoMPositionTolerance = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(0.02);
    settingsStruct.constrainTargetCoMPositionRange.setTimeInterval(settingsStruct.horizon * 0.6, settingsStruct.horizon);

    settingsStruct.comPositionConstraintTolerance = 1e-5*0;
    settingsStruct.centroidalMomentumConstraintTolerance = 1e-5*0;
    settingsStruct.quaternionModulusConstraintTolerance = 1e-2*0;
    settingsStruct.pointPositionConstraintTolerance = 1e-4*0;

    iDynTree::toEigen(settingsStruct.forceMaximumDerivative).setConstant(100.0);
    settingsStruct.normalForceDissipationRatio = 200.0;
//    settingsStruct.normalForceHyperbolicSecantScaling = 300.0;

    //ContactFrictionConstraint
    settingsStruct.frictionCoefficient = 0.3;

    settingsStruct.minimumFeetDistance = 0.10;

    //ContactVelocityControlConstraints
    iDynTree::toEigen(settingsStruct.velocityMaximumDerivative).setConstant(10.0);
    settingsStruct.velocityMaximumDerivative(0) = 10.0;
    settingsStruct.velocityMaximumDerivative(1) = 10.0;
    settingsStruct.planarVelocityHyperbolicTangentScaling = 10.0; //scales the position along z
//    settingsStruct.normalVelocityHyperbolicSecantScaling = 1.0; //scales the force along z

    settingsStruct.complementarityDissipation = 10.0;

    settingsStruct.minimumCoMHeight = 0.5 * initialState.comPosition(2);

    ok = settings.setFromStruct(settingsStruct);
    ASSERT_IS_TRUE(ok);

    auto ipoptSolver = std::make_shared<iDynTree::optimization::IpoptInterface>();

    ASSERT_IS_TRUE(ipoptSolver->isAvailable());

    ok = ipoptSolver->setIpoptOption("linear_solver", "ma97");
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("ma97_order", "metis");
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
    ok = ipoptSolver->setIpoptOption("expect_infeasible_problem", "yes");
    ok = ipoptSolver->setIpoptOption("required_infeasibility_reduction", 0.8);

//    ok = ipoptSolver->setIpoptOption("limited_memory_aug_solver", "extended");
//    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("mu_strategy", "adaptive");
//    ASSERT_IS_TRUE(ok);

//    auto eulerIntegrator = std::make_shared<iDynTree::optimalcontrol::integrators::ForwardEuler>();

//    ok = solver.setIntegrator(eulerIntegrator);
//    ASSERT_IS_TRUE(ok);

    auto worhpSolver = std::make_shared<iDynTree::optimization::WorhpInterface>();

    worhpSolver->setWorhpParam("TolOpti", 1e-4);
    worhpSolver->setWorhpParam("TolFeas", 1e-4);
    worhpSolver->setWorhpParam("TolComp", 1e-5);
    worhpSolver->setWorhpParam("AcceptTolOpti", 1e-1);
    worhpSolver->setWorhpParam("AcceptTolFeas", 1e-3);
//    worhpSolver->setWorhpParam("Algorithm", 1);
//    worhpSolver->setWorhpParam("LineSearchMethod", 3);
//    worhpSolver->setWorhpParam("ArmijoMinAlpha", 5.0e-10);
//    worhpSolver->setWorhpParam("ArmijoMinAlphaRec", 1e-9);

//    worhpSolver->setWorhpParam("Crossover", 2);
//    worhpSolver->setWorhpParam("FeasibleDual", true);



//    ok = solver.setOptimizer(worhpSolver);
//    ASSERT_IS_TRUE(ok);

    ok = solver.setOptimizer(ipoptSolver);
    ASSERT_IS_TRUE(ok);

    ok = solver.specifySettings(settings);
    ASSERT_IS_TRUE(ok);

    ok = solver.setInitialState(initialState);
    ASSERT_IS_TRUE(ok);

    auto stateGuesses = std::make_shared<StateGuess>(comReference, initialState);
    auto controlGuesses = std::make_shared<DynamicalPlanner::TimeInvariantControl>(DynamicalPlanner::Control(vectorList.size(), settingsStruct.leftPointsPosition.size()));

    ok = solver.setGuesses(stateGuesses, controlGuesses);

    std::vector<DynamicalPlanner::State> optimalStates;
    std::vector<DynamicalPlanner::Control> optimalControls;

    DynamicalPlanner::Visualizer visualizer;
    ok = visualizer.setModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
    visualizer.setCameraPosition(iDynTree::Position(1.5, 0.0, 0.0));
    ok = visualizer.visualizeState(initialState);
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
    timeString << timeStruct.tm_year + 1900 << "-" << timeStruct.tm_mon << "-";
    timeString << timeStruct.tm_mday << "_" << timeStruct.tm_hour << "_" << timeStruct.tm_min;
    timeString << "_" << timeStruct.tm_sec;

    ok = visualizer.visualizeStatesAndSaveAnimation(optimalStates, getAbsDirPath("SavedVideos"), "test-1stIteration-" + timeString.str(), "gif", settingsStruct.horizon * settingsStruct.activeControlPercentage);
    ASSERT_IS_TRUE(ok);

//    ok = visualizer.visualizeStates(optimalStates, settingsStruct.horizon * settingsStruct.activeControlPercentage);
//    ASSERT_IS_TRUE(ok);

    //ok = ipoptSolver->setIpoptOption("print_level", 3);
    //ASSERT_IS_TRUE(ok);

//    begin = std::chrono::steady_clock::now();
//    ok = solver.solve(optimalStates, optimalControls);
//    ASSERT_IS_TRUE(ok);
//    end= std::chrono::steady_clock::now();
//    std::cout << "Elapsed time (2nd): " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0 <<std::endl;

//    ok = visualizer.visualizeStates(optimalStates, settingsStruct.horizon * settingsStruct.activeControlPercentage);
//    ASSERT_IS_TRUE(ok);

//    begin = std::chrono::steady_clock::now();
//    ok = solver.solve(optimalStates, optimalControls);
//    ASSERT_IS_TRUE(ok);
//    end= std::chrono::steady_clock::now();
//    std::cout << "Elapsed time (3rd): " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0 <<std::endl;

//    auto timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
//    struct tm timeStruct;
//    timeStruct = *std::localtime(&timeNow);
//    std::ostringstream timeString;
//    timeString << timeStruct.tm_year + 1900 << "-" << timeStruct.tm_mon << "-";
//    timeString << timeStruct.tm_mday << "_" << timeStruct.tm_hour << "_" << timeStruct.tm_min;
//    timeString << "_" << timeStruct.tm_sec;

//    ok = visualizer.visualizeStatesAndSaveAnimation(optimalStates, getAbsDirPath("SavedVideos"), "test-" + timeString.str(), "gif", settingsStruct.horizon * settingsStruct.activeControlPercentage);
//    ASSERT_IS_TRUE(ok);


    ok = ipoptSolver->setIpoptOption("print_level", 0);
    ASSERT_IS_TRUE(ok);

    std::vector<DynamicalPlanner::State> mpcStates;
    mpcStates.push_back(initialState);

    for (size_t i = 0; i < 50; ++i) {
        double initialTime;
        initialState = mpcStates.back();
        initialTime = initialState.time;
        reconstructState(kinDyn, settingsStruct, initialState);
        ok = solver.setInitialState(initialState);
        ASSERT_IS_TRUE(ok);
        iDynTree::toEigen(comReference->get()) += iDynTree::toEigen(iDynTree::Position(0.005, 0.005, 0.0));
        iDynTree::toEigen(meanPointReferencePointer->get()) += iDynTree::toEigen(iDynTree::Position(0.00, 0.00, 0.0));
        begin = std::chrono::steady_clock::now();
        ok = solver.solve(optimalStates, optimalControls);
        if (!ok)
            break;
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time (" << i << "): " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0 <<std::endl;
        optimalStates.front().time += initialTime;
        mpcStates.push_back(optimalStates.front());
        visualizer.visualizeState(mpcStates.back());
    }

    timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    timeStruct = *std::localtime(&timeNow);
    timeString << timeStruct.tm_year + 1900 << "-" << timeStruct.tm_mon << "-";
    timeString << timeStruct.tm_mday << "_" << timeStruct.tm_hour << "_" << timeStruct.tm_min;
    timeString << "_" << timeStruct.tm_sec;

    ok = visualizer.visualizeStatesAndSaveAnimation(mpcStates, getAbsDirPath("SavedVideos"), "test-" + timeString.str(), "gif", settingsStruct.horizon * settingsStruct.activeControlPercentage);
    ASSERT_IS_TRUE(ok);

    return EXIT_SUCCESS;
}
