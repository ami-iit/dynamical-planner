/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#include <iDynTree/Optimizer.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/Utils.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <Eigen/Dense>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlanner/Solver.h>
#include <DynamicalPlanner/RectangularFoot.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <iDynTree/KinDynComputations.h>
#include <URDFdir.h>
#include <chrono>

class OptimizerTest : public iDynTree::optimization::Optimizer {

public:
    OptimizerTest() {}

    virtual ~OptimizerTest() override;

    virtual bool isAvailable() const override{
        return true;
    }

    virtual bool solve() override {
        iDynTree::VectorDynSize variables, constraintsLB, constraintsUB, variablesLB, variablesUB, constraints;
        iDynTree::MatrixDynSize dummyMatrix, jacobian;
        std::vector<size_t> nnzeroRows, nnzeroCols;
        ASSERT_IS_TRUE(m_problem != nullptr);
        ASSERT_IS_TRUE(m_problem->prepare());
        variables.resize(m_problem->numberOfVariables());
//        iDynTree::toEigen(dummyVariables).setConstant(1.0);
        iDynTree::getRandomVector(variables);
        ASSERT_IS_TRUE(m_problem->getConstraintsBounds(constraintsLB, constraintsUB));
        ASSERT_IS_TRUE(constraintsLB.size() == constraintsUB.size());
        for (unsigned int i = 0; i < constraintsLB.size(); ++i) {
            ASSERT_IS_TRUE(constraintsLB(i) <= constraintsUB(i));
        }

        ASSERT_IS_TRUE(m_problem->getVariablesUpperBound(variablesUB));
        ASSERT_IS_TRUE(m_problem->getVariablesLowerBound(variablesLB));
        ASSERT_IS_TRUE(variablesLB.size() == variablesUB.size());
        for (unsigned int i = 0; i < variablesLB.size(); ++i) {
            ASSERT_IS_TRUE(variablesLB(i) <= variablesUB(i));
        }

        ASSERT_IS_TRUE(m_problem->getConstraintsJacobianInfo(nnzeroRows, nnzeroCols));
        ASSERT_IS_TRUE(nnzeroRows.size() == nnzeroCols.size());

        ASSERT_IS_TRUE(m_problem->setVariables(variables));

        jacobian.resize(m_problem->numberOfConstraints(), m_problem->numberOfVariables());
        jacobian.zero();
        ASSERT_IS_TRUE(m_problem->evaluateConstraintsJacobian(jacobian));
        dummyMatrix.resize(jacobian.rows(), jacobian.cols());
        dummyMatrix.zero();

        for (size_t i =0; i < nnzeroRows.size(); ++i){
            unsigned int row = static_cast<unsigned int>(nnzeroRows[i]);
            unsigned int col = static_cast<unsigned int>(nnzeroCols[i]);
            jacobian(row, col) = 0;
        }

        ASSERT_EQUAL_MATRIX_TOL(dummyMatrix, jacobian, iDynTree::DEFAULT_TOL); //check the sparsity structure

        ASSERT_IS_TRUE(m_problem->getGuess(variables));

        for (unsigned int i = 0; i < variables.size(); ++i) {
            if (iDynTree::checkDoublesAreEqual(variablesLB(i), variablesUB(i))) {
                ASSERT_EQUAL_DOUBLE_TOL(variables(i), variablesLB(i), iDynTree::DEFAULT_TOL);
            } else {
                ASSERT_IS_TRUE(variables(i) >= variablesLB(i));
                ASSERT_IS_TRUE(variables(i) <= variablesUB(i));
            }
        }

        ASSERT_IS_TRUE(m_problem->setVariables(variables));
        ASSERT_IS_TRUE(m_problem->evaluateConstraints(constraints));

        for (unsigned int i = 0; i < constraints.size(); ++i) {
            if (iDynTree::checkDoublesAreEqual(constraintsLB(i), constraintsUB(i))) {
                ASSERT_EQUAL_DOUBLE_TOL(constraints(i), constraintsLB(i), 1e-3);
            } else {
                ASSERT_IS_TRUE(constraints(i) >= (constraintsLB(i) - 1e-3));
                ASSERT_IS_TRUE(constraints(i) <= (constraintsUB(i) + 1e-3));
            }
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

void fillInitialState(const iDynTree::Model& model, const DynamicalPlanner::SettingsStruct settings,
                      const iDynTree::VectorDynSize desiredJoints, DynamicalPlanner::RectangularFoot &foot,
                      DynamicalPlanner::State &initialState) {

    iDynTree::KinDynComputations kinDyn;

    bool ok = kinDyn.loadRobotModel(model);
    ASSERT_IS_TRUE(ok);

    ok = kinDyn.setFloatingBase(model.getLinkName(model.getFrameLink(model.getFrameIndex(settings.leftFrameName))));
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

    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

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

    iDynTree::toEigen(desiredInitialJoints) *= iDynTree::deg2rad(1.0);

    fillInitialState(modelLoader.model(), settingsStruct, desiredInitialJoints, foot, initialState);

    auto comReference = std::make_shared<CoMReference>(initialState.comPosition, 0.0, 0.0);

    settingsStruct.desiredCoMTrajectory  = comReference;

    settingsStruct.desiredJointsTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredInitialJoints);

    settingsStruct.minimumDt = 0.1;
    settingsStruct.controlPeriod = 0.1;
    settingsStruct.maximumDt = 1.0;
    settingsStruct.horizon = 1.0;

    ok = settings.setFromStruct(settingsStruct);
    ASSERT_IS_TRUE(ok);

    auto optimizerTest = std::make_shared<OptimizerTest>();


    ok = solver.setOptimizer(optimizerTest);
    ASSERT_IS_TRUE(ok);

    ok = solver.specifySettings(settings);
    ASSERT_IS_TRUE(ok);

    ok = solver.setInitialState(initialState);
    ASSERT_IS_TRUE(ok);

    auto stateGuesses = std::make_shared<StateGuess>(comReference, initialState);
    auto controlGuesses = std::make_shared<DynamicalPlanner::TimeInvariantControl>(DynamicalPlanner::Control(vectorList.size(), settingsStruct.leftPointsPosition.size()));

    std::vector<DynamicalPlanner::State> optimalStates;
    std::vector<DynamicalPlanner::Control> optimalControls;

    for (size_t i = 0; i < 5; ++i) {
        ok = solver.setGuesses(stateGuesses, controlGuesses);
        ok = solver.solve(optimalStates, optimalControls);
        ASSERT_IS_TRUE(ok);
    }

    return EXIT_SUCCESS;
}
