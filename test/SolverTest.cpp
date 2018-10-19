/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Solver.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/KinDynComputations.h>
#include <iDynTree/Optimizers/IpoptInterface.h>
#include <iDynTree/Optimizers/WorhpInterface.h>
#include <iDynTree/Integrators/ForwardEuler.h>
#include <URDFdir.h>

void fillInitialState(const iDynTree::Model& model, const DynamicalPlanner::SettingsStruct settings,
                      const iDynTree::VectorDynSize desiredJoints, DynamicalPlanner::State &initialState) {

    iDynTree::KinDynComputations kinDyn;

    bool ok = kinDyn.loadRobotModel(model);
    ASSERT_IS_TRUE(ok);

    ok = kinDyn.setFloatingBase(model.getLinkName(model.getFrameLink(model.getFrameIndex(settings.leftFrameName))));
    ASSERT_IS_TRUE(ok);

    iDynTree::Vector3 gravity;
    gravity.zero();
    gravity(2) = -9.81;

    ok = kinDyn.setRobotState(model.getFrameTransform(model.getFrameIndex(settings.leftFrameName)).inverse(), desiredJoints,
                              iDynTree::Twist::Zero(), iDynTree::VectorDynSize(desiredJoints.size()), gravity);
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

    double pointForce = totalMass * 9.81 / (2 * settings.leftPointsPosition.size());

    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointPosition = leftTransform * settings.leftPointsPosition[i];
        initialState.rightContactPointsState[i].pointPosition = rightTransform * settings.rightPointsPosition[i];

        initialState.leftContactPointsState[i].pointForce.zero();
        initialState.leftContactPointsState[i].pointForce(2) = pointForce;

        initialState.rightContactPointsState[i].pointForce.zero();
        initialState.rightContactPointsState[i].pointForce(2) = pointForce;
    }
}

class CoMReference : public iDynTree::optimalcontrol::TimeVaryingVector {
    iDynTree::VectorDynSize desiredCoM;
    double xVelocity, yVelocity;
    iDynTree::Vector3 initialCoM;

public:
    CoMReference(iDynTree::Vector3 &CoMinitial, double velX, double velY)
        : desiredCoM(3)
        , xVelocity(velX)
        , yVelocity(velY)
        , initialCoM(CoMinitial)
    { }

    ~CoMReference() override;

    iDynTree::VectorDynSize &get(double time, bool &isValid) override {
        desiredCoM(0) = initialCoM(0) + xVelocity * time;
        desiredCoM(1) = initialCoM(1) + yVelocity * time;
        desiredCoM(2) = initialCoM(2);
        isValid = true;
        return desiredCoM;
    }

};
CoMReference::~CoMReference() { }

class StateGuess : public DynamicalPlanner::TimeVaryingState {
    DynamicalPlanner::State m_state, m_initialState;
    std::shared_ptr<CoMReference> m_comReference;
public:

    StateGuess(std::shared_ptr<CoMReference> comReference, const DynamicalPlanner::State &initialState)
        : m_state(initialState)
        , m_initialState(initialState)
        , m_comReference(comReference)
    { }

    ~StateGuess() override;

    DynamicalPlanner::State &get(double time, bool &isValid) override {
        m_state.comPosition = m_comReference->get(time, isValid);
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

//    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
//                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
//                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
//                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "r_shoulder_pitch", "r_shoulder_roll", "l_hip_pitch", "l_hip_roll", "l_hip_yaw",
                                         "l_knee", "l_ankle_pitch", "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw",
                                         "r_knee", "r_ankle_pitch", "r_ankle_roll"});

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


    DynamicalPlanner::State initialState;

    initialState.resize(vectorList.size(), settingsStruct.leftPointsPosition.size());

    iDynTree::VectorDynSize desiredInitialJoints(static_cast<unsigned int>(modelLoader.model().getNrOfDOFs()));

//    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, -7, 22, 11, 30, -7, 22, 11, 30, 5.082, 0.406, -0.131,
//                                              -45.249, -26.454, -0.351, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351;

    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, -7, 22, -7, 22, 5.082, 0.406, -0.131, -45.249,
                                             -26.454, -0.351, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351;

//    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351, 5.082,
//                                               0.406, -0.131, -45.249, -26.454, -0.351;

    iDynTree::toEigen(desiredInitialJoints) *= iDynTree::deg2rad(1.0);


    fillInitialState(modelLoader.model(), settingsStruct, desiredInitialJoints, initialState);

    auto comReference = std::make_shared<CoMReference>(initialState.comPosition, 1.0, 0.0);

    settingsStruct.desiredCoMTrajectory  = comReference;

    settingsStruct.desiredJointsTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredInitialJoints);

    settingsStruct.frameCostActive = true;
    settingsStruct.staticTorquesCostActive = true;
    settingsStruct.forceMeanCostActive = true;
    settingsStruct.comCostActive = true;
    settingsStruct.forceDerivativeCostActive = true;
    settingsStruct.pointAccelerationCostActive = true;
    settingsStruct.jointsRegularizationCostActive = true;
    settingsStruct.jointsVelocityCostActive = true;

    settingsStruct.jointsVelocityCostOverallWeight = 1e-4;
    settingsStruct.staticTorquesCostOverallWeight = 1e-5;

//    settingsStruct.minimumDt = 0.01;
//    settingsStruct.controlPeriod = 0.1;
//    settingsStruct.maximumDt = 0.05;

    settingsStruct.minimumDt = 0.1;
    settingsStruct.controlPeriod = 0.1;
    settingsStruct.maximumDt = 1.0;

    settingsStruct.comPositionConstraintTolerance = 1e-2;
    settingsStruct.centroidalMomentumConstraintTolerance = 1e-2;
    settingsStruct.quaternionModulusConstraintTolerance = 1e-2;
    settingsStruct.pointPositionConstraintTolerance = 5e-3;

    iDynTree::toEigen(settingsStruct.forceMaximumDerivative).setConstant(100.0);
    settingsStruct.normalForceDissipationRatio = 1.0;
    settingsStruct.normalForceHyperbolicSecantScaling = 200.0;

    //ContactFrictionConstraint
    settingsStruct.frictionCoefficient = 0.3;

    settingsStruct.minimumFeetDistance = 0.05;

    //ContactVelocityControlConstraints
    iDynTree::toEigen(settingsStruct.velocityMaximumDerivative).setConstant(10.0);
    settingsStruct.planarVelocityHyperbolicTangentScaling = 10.0; //scales the position along z
    settingsStruct.normalVelocityHyperbolicSecantScaling = 1.0; //scales the force along z

    ok = settings.setFromStruct(settingsStruct);
    ASSERT_IS_TRUE(ok);

    auto ipoptSolver = std::make_shared<iDynTree::optimization::IpoptInterface>();

    ASSERT_IS_TRUE(ipoptSolver->isAvailable());

    ok = ipoptSolver->setIpoptOption("linear_solver", "ma27");
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("print_level", 5);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("nlp_scaling_max_gradient", 10.0);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("nlp_scaling_min_value", 1e-6);
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("evaluate_orig_obj_at_resto_trial", "no");
//    ASSERT_IS_TRUE(ok);

//    ok = ipoptSolver->setIpoptOption("check_derivatives_for_naninf", "yes");
//    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("tol", 1e-5);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("constr_viol_tol", 1e-3);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("acceptable_tol", 1e-3);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("acceptable_iter", 3);
    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("alpha_for_y", "min-dual-infeas");
    ASSERT_IS_TRUE(ok);
//    ok = ipoptSolver->setIpoptOption("mu_strategy", "adaptive");
//    ASSERT_IS_TRUE(ok);

//    auto eulerIntegrator = std::make_shared<iDynTree::optimalcontrol::integrators::ForwardEuler>();

//    ok = solver.setIntegrator(eulerIntegrator);
//    ASSERT_IS_TRUE(ok);

    auto worhpSolver = std::make_shared<iDynTree::optimization::WorhpInterface>();

    worhpSolver->setWorhpParam("TolOpti", 1e-4);
    worhpSolver->setWorhpParam("TolFeas", 1e-4);
    worhpSolver->setWorhpParam("TolComp", 1e-5);
    worhpSolver->setWorhpParam("AcceptTolOpti", 1e-3);
    worhpSolver->setWorhpParam("AcceptTolFeas", 1e-3);
    worhpSolver->setWorhpParam("Algorithm", 1);
    worhpSolver->setWorhpParam("LineSearchMethod", 3);
    worhpSolver->setWorhpParam("ArmijoMinAlpha", 5.0e-10);
    worhpSolver->setWorhpParam("ArmijoMinAlphaRec", 1e-9);

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

    ok = solver.solve(optimalStates, optimalControls);

    ASSERT_IS_TRUE(ok);
    ok = ipoptSolver->setIpoptOption("print_level", 3);
    ASSERT_IS_TRUE(ok);

    ok = solver.solve(optimalStates, optimalControls);
    ASSERT_IS_TRUE(ok);


    return EXIT_SUCCESS;
}