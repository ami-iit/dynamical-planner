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
#include <DynamicalPlanner/Utilities.h>
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

    DynamicalPlanner::Utilities::FillDefaultInitialState(settingsStruct, desiredInitialJoints, foot, foot, initialState);

    iDynTree::VectorDynSize comPointReference(3);
    comPointReference = initialState.comPosition;
    auto comReference = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comPointReference);

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

    auto stateGuesses = std::make_shared<DynamicalPlanner::Utilities::TranslatingCoMStateGuess>(comReference, initialState);
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
