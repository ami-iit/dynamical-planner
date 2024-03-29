/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Solver.h>
#include <DynamicalPlanner/RectangularFoot.h>
#include <DynamicalPlanner/Utilities.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/Integrators/ForwardEuler.h>
#include <URDFdir.h>
#include <FolderPath.h>
#include <chrono>
#include <ctime>
#include <sstream>

class OptimizerTest : public iDynTree::optimization::Optimizer {

public:
    OptimizerTest() {}

    virtual ~OptimizerTest() override;

    virtual bool isAvailable() const override{
        return true;
    }

    virtual bool solve() override {
        iDynTree::VectorDynSize variables, constraints, costGradient;
        double cost;
        iDynTree::MatrixDynSize jacobian;
        ASSERT_IS_TRUE(m_problem != nullptr);
        ASSERT_IS_TRUE(m_problem->prepare());
        variables.resize(m_problem->numberOfVariables());
//        iDynTree::toEigen(dummyVariables).setConstant(1.0);
        jacobian.resize(m_problem->numberOfConstraints(), m_problem->numberOfVariables());
        jacobian.zero();
        constraints.resize(m_problem->numberOfConstraints());
        costGradient.resize(m_problem->numberOfVariables());

        for (size_t i = 0; i < 1000 ; ++i) {
            iDynTree::getRandomVector(variables);
            ASSERT_IS_TRUE(m_problem->setVariables(variables));
            ASSERT_IS_TRUE(m_problem->evaluateConstraints(constraints));
            ASSERT_IS_TRUE(m_problem->evaluateConstraintsJacobian(jacobian));
            ASSERT_IS_TRUE(m_problem->evaluateCostFunction(cost));
            ASSERT_IS_TRUE(m_problem->evaluateCostGradient(costGradient));
        }

        return true;
    }

    virtual bool getPrimalVariables(iDynTree::VectorDynSize &primalVariables) override {
        ASSERT_IS_TRUE(m_problem != nullptr);
        primalVariables.resize(m_problem->numberOfVariables());
        primalVariables.zero();
        return true;
    }

    virtual bool getDualVariables(iDynTree::VectorDynSize &constraintsMultipliers,
                                  iDynTree::VectorDynSize &lowerBoundsMultipliers,
                                  iDynTree::VectorDynSize &upperBoundsMultipliers) override {
        ASSERT_IS_TRUE(m_problem != nullptr);
        constraintsMultipliers.resize(m_problem->numberOfConstraints());
        lowerBoundsMultipliers.resize(m_problem->numberOfVariables());
        upperBoundsMultipliers.resize(m_problem->numberOfVariables());
        return true;
    }
};
OptimizerTest::~OptimizerTest(){}

int main() {

    DynamicalPlanner::Solver solver;

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

//    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, -7, 22, 11, 30, -7, 22, 11, 30, 5.082, 0.406, -0.131,
//                                              -45.249, -26.454, -0.351, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351;

    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, -7, 22, -7, 22, 5.082, 0.406, -0.131, -45.249,
                                             -26.454, -0.351, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351;

//    iDynTree::toEigen(desiredInitialJoints) << 15, 0, 0, 5.082, 0.406, -0.131, -45.249, -26.454, -0.351, 5.082,
//                                               0.406, -0.131, -45.249, -26.454, -0.351;

    iDynTree::toEigen(desiredInitialJoints) *= iDynTree::deg2rad(1.0);


    DynamicalPlanner::Utilities::FillDefaultInitialState(settingsStruct, desiredInitialJoints, foot, foot, initialState);

//    auto comReference = std::make_shared<CoMReference>(initialState.comPosition, 0.2, 0.0, 0.0);
    iDynTree::VectorDynSize comPointReference(3);
    iDynTree::toEigen(comPointReference) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.3, 0.1, 0.0));
    auto comReference = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comPointReference);

    settingsStruct.desiredCoMTrajectory  = comReference;

    settingsStruct.desiredJointsTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredInitialJoints);

    settingsStruct.frameCostActive = true;
    settingsStruct.staticTorquesCostActive = false;
    settingsStruct.forceMeanCostActive = false;
    settingsStruct.comCostActive = false;
    settingsStruct.forceDerivativeCostActive = true;
    settingsStruct.pointAccelerationCostActive = true;
    settingsStruct.jointsRegularizationCostActive = true;
    settingsStruct.jointsVelocityCostActive = true;
    settingsStruct.swingCostActive = true;
    settingsStruct.phantomForcesCostActive = false;

    settingsStruct.jointsVelocityCostOverallWeight = 1e-4;
    settingsStruct.staticTorquesCostOverallWeight = 1e-5;
    settingsStruct.jointsRegularizationCostOverallWeight = 1e-1;
    settingsStruct.jointsVelocityCostOverallWeight = 1e-20;
    settingsStruct.forceMeanCostOverallWeight = 1.0;
    settingsStruct.forceDerivativesCostOverallWeight = 1e-15;
    settingsStruct.pointAccelerationCostOverallWeight = 1e-15;
    settingsStruct.swingCostOverallWeight = 10;
    settingsStruct.phantomForcesCostOverallWeight = 1.0;
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
    settingsStruct.horizon = 2.0;
    settingsStruct.activeControlPercentage = 1.0;
    settingsStruct.comCostActiveRange.setTimeInterval(settingsStruct.horizon * settingsStruct.activeControlPercentage, settingsStruct.horizon);

    settingsStruct.constrainTargetCoMPosition = true;
    settingsStruct.targetCoMPositionTolerance = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(0.02);
    settingsStruct.constrainTargetCoMPositionRange.setTimeInterval(settingsStruct.horizon * 0.8, settingsStruct.horizon);

    settingsStruct.comPositionConstraintTolerance = 1e-5;
    settingsStruct.centroidalMomentumConstraintTolerance = 1e-5;
    settingsStruct.quaternionModulusConstraintTolerance = 1e-2;
    settingsStruct.pointPositionConstraintTolerance = 1e-4;

    iDynTree::toEigen(settingsStruct.forceMaximumDerivative).setConstant(50.0);
    settingsStruct.normalForceDissipationRatio = 50.0;
    settingsStruct.normalForceHyperbolicSecantScaling = 150.0;

    //ContactFrictionConstraint
    settingsStruct.frictionCoefficient = 0.3;

    settingsStruct.minimumFeetDistance = 0.05;

    //ContactVelocityControlConstraints
    iDynTree::toEigen(settingsStruct.velocityMaximumDerivative).setConstant(10.0);
    settingsStruct.velocityMaximumDerivative(0) = 10.0;
    settingsStruct.velocityMaximumDerivative(1) = 10.0;
    settingsStruct.planarVelocityHyperbolicTangentScaling = 10.0; //scales the position along z
    settingsStruct.normalVelocityHyperbolicSecantScaling = 1.0; //scales the force along z

    settingsStruct.complementarityDissipation = 10.0;

    settingsStruct.minimumCoMHeight = 0.95 * initialState.comPosition(2);

    auto optimizerTest = std::make_shared<OptimizerTest>();

    ok = solver.setOptimizer(optimizerTest);
    ASSERT_IS_TRUE(ok);

    ok = solver.specifySettings(settingsStruct);
    ASSERT_IS_TRUE(ok);

    ok = solver.setInitialState(initialState);
    ASSERT_IS_TRUE(ok);

    auto stateGuesses = std::make_shared<DynamicalPlanner::Utilities::TranslatingCoMStateGuess>(comReference, initialState);
    auto controlGuesses = std::make_shared<DynamicalPlanner::TimeInvariantControl>(DynamicalPlanner::Control(vectorList.size(), settingsStruct.leftPointsPosition.size()));

    ok = solver.setGuesses(stateGuesses, controlGuesses);

    std::vector<DynamicalPlanner::State> optimalStates;
    std::vector<DynamicalPlanner::Control> optimalControls;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    ok = solver.solve(optimalStates, optimalControls);
    ASSERT_IS_TRUE(ok);
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time (1st): " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count())/1000.0 <<std::endl;

    return EXIT_SUCCESS;
}
