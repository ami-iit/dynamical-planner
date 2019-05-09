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
#include <iDynTree/Integrators/ForwardEuler.h>
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
        iDynTree::VectorDynSize variables, constraints, multipliers;
        iDynTree::MatrixDynSize dummyMatrix, jacobian, constraintsHessian, costHessian, hessian;
        std::vector<size_t> nnzeroRowsJac, nnzeroColsJac, nnZeroRowsHes, nnZeroColsHes;
        ASSERT_IS_TRUE(m_problem != nullptr);
        ASSERT_IS_TRUE(m_problem->prepare());
        variables.resize(m_problem->numberOfVariables());
        iDynTree::getRandomVector(variables);

        ASSERT_IS_TRUE(m_problem->setVariables(variables));

        ASSERT_IS_TRUE(m_problem->getConstraintsJacobianInfo(nnzeroRowsJac, nnzeroColsJac));
        ASSERT_IS_TRUE(nnzeroRowsJac.size() == nnzeroColsJac.size());

        jacobian.resize(m_problem->numberOfConstraints(), m_problem->numberOfVariables());
        jacobian.zero();
        ASSERT_IS_TRUE(m_problem->evaluateConstraintsJacobian(jacobian));
        dummyMatrix.resize(jacobian.rows(), jacobian.cols());
        dummyMatrix.zero();

        for (size_t i =0; i < nnzeroRowsJac.size(); ++i){
            unsigned int row = static_cast<unsigned int>(nnzeroRowsJac[i]);
            unsigned int col = static_cast<unsigned int>(nnzeroColsJac[i]);
            jacobian(row, col) = 0;
        }

        ASSERT_EQUAL_MATRIX_TOL(dummyMatrix, jacobian, iDynTree::DEFAULT_TOL); //check the sparsity structure

        multipliers.resize(m_problem->numberOfConstraints());

        iDynTree::getRandomVector(multipliers);

        costHessian.resize(m_problem->numberOfVariables(), m_problem->numberOfVariables());
        costHessian.zero();

        constraintsHessian.resize(m_problem->numberOfVariables(), m_problem->numberOfVariables());
        constraintsHessian.zero();

        ASSERT_IS_TRUE(m_problem->evaluateCostHessian(costHessian));

        ASSERT_IS_TRUE(m_problem->evaluateConstraintsHessian(multipliers, constraintsHessian));

        auto costHessianEigen = iDynTree::toEigen(costHessian);
        auto constraintsHessianEigen = iDynTree::toEigen(constraintsHessian);

        hessian.resize(costHessian.rows(), costHessian.cols());

        iDynTree::toEigen(hessian) = costHessianEigen + constraintsHessianEigen;

        ASSERT_IS_TRUE(m_problem->getHessianInfo(nnZeroRowsHes, nnZeroColsHes));

        dummyMatrix.resize(hessian.rows(), hessian.cols());
        dummyMatrix.zero();

        for (size_t i =0; i < nnZeroRowsHes.size(); ++i){
            unsigned int row = static_cast<unsigned int>(nnZeroRowsHes[i]);
            unsigned int col = static_cast<unsigned int>(nnZeroColsHes[i]);
            hessian(row, col) = 0;
        }

        ASSERT_EQUAL_MATRIX_TOL(dummyMatrix, hessian, iDynTree::DEFAULT_TOL); //check the sparsity structure

        bool ok = m_problem->setVariables(variables);
        ASSERT_IS_TRUE(ok);

        ASSERT_IS_TRUE(m_problem->evaluateCostHessian(costHessian));

        ASSERT_IS_TRUE(m_problem->evaluateConstraintsHessian(multipliers, constraintsHessian));

        iDynTree::toEigen(hessian) = costHessianEigen + constraintsHessianEigen;

        for (unsigned int i = 0; i < hessian.rows(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                ASSERT_EQUAL_DOUBLE(hessian(i, j), hessian(j, i)); //check symmetricity
            }
        }


        //-------------------------------- Check Constraints hessian
        {
            double perturbation = 1e-4;
            iDynTree::VectorDynSize firstOrderTaylor, perturbedRow, originalVariables, perturbedVariables, perturbedConstraints;
            iDynTree::MatrixDynSize originalJacobian = jacobian;
            originalVariables = variables;


            iDynTree::MatrixDynSize perturbedJacobian;
            perturbedJacobian = jacobian;

            ok = m_problem->setVariables(originalVariables);
            ASSERT_IS_TRUE(ok);

            ok = m_problem->evaluateConstraintsJacobian(originalJacobian);
            ASSERT_IS_TRUE(ok);

            iDynTree::VectorDynSize originalConstraints(multipliers.size());
            ok = m_problem->evaluateConstraints(originalConstraints);

            perturbedConstraints = originalConstraints;
            perturbedVariables = originalVariables;

            for (unsigned int row = 0; row < multipliers.size(); ++row) {
                multipliers.zero();
                multipliers(row) = 1.0;

                std::cerr << "Row: " << row << std::endl;

                constraintsHessian.zero();

                perturbedVariables = originalVariables;
                ok = m_problem->setVariables(originalVariables);
                ASSERT_IS_TRUE(ok);

                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                ok = m_problem->evaluateConstraintsHessian(multipliers, constraintsHessian);
                ASSERT_IS_TRUE(ok);
                std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
                if (!row) {
                    std::cout << "Elapsed time ms hessian: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
                }

                perturbedRow.resize(originalJacobian.cols());
                firstOrderTaylor = perturbedRow;

                for (unsigned int i = 0; i < originalVariables.size(); ++i) {
                    //                std::cerr << "Variable: " << i << std::endl;

                    perturbedVariables = originalVariables;
                    perturbedVariables(i) = perturbedVariables(i) + perturbation;

                    ok = m_problem->setVariables(perturbedVariables);
                    ASSERT_IS_TRUE(ok);

                    ok = m_problem->evaluateConstraintsJacobian(perturbedJacobian);
                    ASSERT_IS_TRUE(ok);

                    iDynTree::toEigen(perturbedRow) = iDynTree::toEigen(perturbedJacobian).row(row).transpose();

                    iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalJacobian).row(row).transpose() +
                        iDynTree::toEigen(constraintsHessian) * (iDynTree::toEigen(perturbedVariables) - iDynTree::toEigen(originalVariables));
                    ASSERT_EQUAL_VECTOR_TOL(perturbedRow, firstOrderTaylor, perturbation/10);
                }
            }

            for (size_t i = 0; i < 3; ++i) {
                ok = m_problem->setVariables(originalVariables);
                ASSERT_IS_TRUE(ok);

                iDynTree::getRandomVector(multipliers,-10.0, 10.0);

                ok = m_problem->evaluateConstraintsHessian(multipliers, constraintsHessian);
                ASSERT_IS_TRUE(ok);

                iDynTree::VectorDynSize delta = originalVariables;

                iDynTree::getRandomVector(delta);

                iDynTree::toEigen(delta).normalize();

                //            iDynTree::toEigen(delta) *= perturbation*100;

                iDynTree::toEigen(perturbedVariables) = iDynTree::toEigen(originalVariables) + iDynTree::toEigen(delta);

                double originalLagrangian = iDynTree::toEigen(multipliers).transpose() * iDynTree::toEigen(originalConstraints);

//                Eigen::JacobiSVD<Eigen::MatrixXd> svd(iDynTree::toEigen(constraintsHessian), Eigen::ComputeThinU | Eigen::ComputeThinV);

//                iDynTree::VectorDynSize deltaProjected = delta;

//                iDynTree::toEigen(deltaProjected) = svd.matrixV() * iDynTree::toEigen(delta);

                double firstOrderApprox = originalLagrangian + (iDynTree::toEigen(delta).transpose() * iDynTree::toEigen(originalJacobian).transpose() * iDynTree::toEigen(multipliers))(0,0);
                //            double secondOrderApprox = firstOrderApprox + 0.5 * iDynTree::toEigen(deltaProjected).topRows(svd.rank()).transpose() * svd.singularValues().topRows(svd.rank()).asDiagonal() * iDynTree::toEigen(deltaProjected).topRows(svd.rank());

                double secondOrderApprox = firstOrderApprox + 0.5 * iDynTree::toEigen(delta).transpose() *
                        iDynTree::toEigen(constraintsHessian) * iDynTree::toEigen(delta);

                iDynTree::MatrixDynSize test = constraintsHessian;

//                iDynTree::toEigen(test) = svd.matrixU() * svd.matrixV();

//                iDynTree::MatrixDynSize testIdentity = constraintsHessian;
//                iDynTree::toEigen(testIdentity).setIdentity();

//                Eigen::MatrixXd thisShouldBeNull = iDynTree::toEigen(test) - iDynTree::toEigen(testIdentity);
//                Eigen::JacobiSVD<Eigen::MatrixXd> svd2(thisShouldBeNull);
//                svd2.setThreshold(1e-10);

//                std::cerr << "Test max eigenvalue: " << svd2.singularValues().maxCoeff() << std::endl;


//                std::cerr << "Hessian dim: " << hessian.rows() << " rank: " << svd.rank() << std::endl;

                std::cerr << "Hessian contribution: " << secondOrderApprox - firstOrderApprox << std::endl;

                ok = m_problem->setVariables(perturbedVariables);
                ASSERT_IS_TRUE(ok);

                ok = m_problem->evaluateConstraints(perturbedConstraints);
                ASSERT_IS_TRUE(ok);

                double perturbedLagrangian = iDynTree::toEigen(multipliers).transpose() * iDynTree::toEigen(perturbedConstraints);

                double firstOrderError = std::abs(perturbedLagrangian - firstOrderApprox);
                std::cerr << "Constraints First order error: " << firstOrderError << std::endl;
                double secondOrderError = std::abs(perturbedLagrangian - secondOrderApprox);
                std::cerr << "Constraint Second order error: " << secondOrderError << std::endl;

                ASSERT_IS_TRUE(firstOrderError >= secondOrderError);
            }

        }

        //-------------------------------- Check Costs hessian
        {
            double perturbation = 1e-4;
            iDynTree::VectorDynSize firstOrderTaylor, perturbedGradient, originalVariables, perturbedVariables;
            iDynTree::VectorDynSize originalGradient;
            iDynTree::MatrixDynSize costHessian;
            originalVariables = variables;

            ok = m_problem->setVariables(originalVariables);
            ASSERT_IS_TRUE(ok);

            originalGradient.resize(originalVariables.size());

            ok = m_problem->evaluateCostGradient(originalGradient);
            ASSERT_IS_TRUE(ok);
            perturbedGradient = originalGradient;

            costHessian.resize(originalVariables.size(), originalVariables.size());

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            ok = m_problem->evaluateCostHessian(costHessian);
            ASSERT_IS_TRUE(ok);
            std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
            std::cout << "Elapsed time ms cost hessian: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
            firstOrderTaylor = originalGradient;


            for (unsigned int i = 0; i < originalVariables.size(); ++i) {
                //                std::cerr << "Variable: " << i << std::endl;

                perturbedVariables = originalVariables;
                perturbedVariables(i) = perturbedVariables(i) + perturbation;

                ok = m_problem->setVariables(perturbedVariables);
                ASSERT_IS_TRUE(ok);

                ok = m_problem->evaluateCostGradient(perturbedGradient);
                ASSERT_IS_TRUE(ok);

                iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalGradient) +
                    iDynTree::toEigen(costHessian) * (iDynTree::toEigen(perturbedVariables) - iDynTree::toEigen(originalVariables));
                ASSERT_EQUAL_VECTOR_TOL(perturbedGradient, firstOrderTaylor, perturbation/10);
            }

            ok = m_problem->setVariables(originalVariables);
            ASSERT_IS_TRUE(ok);

            double originalCost;
            ok = m_problem->evaluateCostFunction(originalCost);

            iDynTree::VectorDynSize delta = originalVariables;

            iDynTree::getRandomVector(delta);

            iDynTree::toEigen(delta).normalize();

            iDynTree::toEigen(perturbedVariables) = iDynTree::toEigen(originalVariables) + iDynTree::toEigen(delta);

            double firstOrderApprox = originalCost + (iDynTree::toEigen(originalGradient).transpose()*iDynTree::toEigen(delta))(0,0);
            double secondOrderApprox = firstOrderApprox + 0.5 * iDynTree::toEigen(delta).transpose() * iDynTree::toEigen(costHessian) * iDynTree::toEigen(delta);
            double perturbedCost;

            ok = m_problem->setVariables(perturbedVariables);
            ASSERT_IS_TRUE(ok);

            ok = m_problem->evaluateCostFunction(perturbedCost);
            ASSERT_IS_TRUE(ok);

            double firstOrderError = std::abs(perturbedCost - firstOrderApprox);
            std::cerr << "Cost First order error: " << firstOrderError << std::endl;
            double secondOrderError = std::abs(perturbedCost - secondOrderApprox);
            std::cerr << "Cost Second order error: " << secondOrderError << std::endl;


            ASSERT_IS_TRUE(firstOrderError >= secondOrderError);

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

typedef struct {
    iDynTree::Position desiredPosition;
    iDynTree::optimalcontrol::TimeRange activeRange;
} PositionWithTimeRange;

class MeanPointReferenceGenerator;

class MeanPointReferenceGeneratorData {

    friend class MeanPointReferenceGenerator;

    MeanPointReferenceGeneratorData(size_t desiredPoints) {

        desiredPositions.resize(desiredPoints);

        for (auto& el : desiredPositions) {
            el.desiredPosition.zero();
        }

    }

public:

    ~MeanPointReferenceGeneratorData() {}

    std::vector<PositionWithTimeRange> desiredPositions;

};

class TimeVaryingWeight : public iDynTree::optimalcontrol::TimeVaryingDouble {

    friend class MeanPointReferenceGenerator;

    std::shared_ptr<MeanPointReferenceGeneratorData> m_data;
    double m_increaseFactor;
    double m_outputWeight;

    TimeVaryingWeight(std::shared_ptr<MeanPointReferenceGeneratorData> data, double increaseFactor) {
        m_data = data;
        m_increaseFactor = increaseFactor;
    }

public:

    ~TimeVaryingWeight() override;


    const double& get(double time, bool& isValid) override {

        isValid = true;
        std::vector<PositionWithTimeRange>::reverse_iterator activeElement;
        activeElement = std::find_if(m_data->desiredPositions.rbegin(),
                                     m_data->desiredPositions.rend(),
                                     [time](const PositionWithTimeRange& a) -> bool { return a.activeRange.isInRange(time); }); //find the last element in the vector with init time lower than the specified time

        if (activeElement == m_data->desiredPositions.rend()) {
            m_outputWeight = 0.0;
            return m_outputWeight;
        }

        double timeWeight = m_increaseFactor * (time - activeElement->activeRange.initTime())/(activeElement->activeRange.endTime() - activeElement->activeRange.initTime()) + 1.0;
        m_outputWeight = timeWeight*timeWeight;
        return m_outputWeight;
    }
};
TimeVaryingWeight::~TimeVaryingWeight(){}

class MeanPointReferencePosition : public iDynTree::optimalcontrol::TimeVaryingPosition {

    friend class MeanPointReferenceGenerator;
    std::shared_ptr<MeanPointReferenceGeneratorData> m_data;
    iDynTree::Position m_zeroPosition;

    MeanPointReferencePosition(std::shared_ptr<MeanPointReferenceGeneratorData> data) {
        m_data = data;
    }

public:

    ~MeanPointReferencePosition() override;

    const iDynTree::Position& get(double time, bool& isValid) override {

        isValid = true;
        if (!m_data->desiredPositions.size()) {
            return m_zeroPosition;
        }

        std::vector<PositionWithTimeRange>::reverse_iterator activeElement;
        activeElement = std::find_if(m_data->desiredPositions.rbegin(),
                                     m_data->desiredPositions.rend(),
                                     [time](const PositionWithTimeRange& a) -> bool { return a.activeRange.isInRange(time); }); //find the last element in the vector with init time lower than the specified time

        if (activeElement == m_data->desiredPositions.rend()) {
            return m_zeroPosition;
        }

        return activeElement->desiredPosition;
    }
};
MeanPointReferencePosition::~MeanPointReferencePosition(){}

class MeanPointReferenceGenerator {

    std::shared_ptr<TimeVaryingWeight> m_weightPointer;
    std::shared_ptr<MeanPointReferencePosition> m_positionPointer;
    std::shared_ptr<MeanPointReferenceGeneratorData> m_data;

public:

    MeanPointReferenceGenerator(unsigned int desiredPoints, double increaseFactor) {
        m_data.reset(new MeanPointReferenceGeneratorData(desiredPoints));
        m_weightPointer.reset(new TimeVaryingWeight(m_data, increaseFactor));
        m_positionPointer.reset(new MeanPointReferencePosition(m_data));
    }

    ~MeanPointReferenceGenerator() {}

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> timeVaryingWeight() {
        return m_weightPointer;
    }

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> timeVaryingReference() {
        return m_positionPointer;
    }

    PositionWithTimeRange& operator[](size_t index) {
        return m_data->desiredPositions[index];
    }

    const PositionWithTimeRange& operator[](size_t index) const {
        return m_data->desiredPositions[index];
    }


    void resize(size_t newSize) {
        PositionWithTimeRange zeroElement;
        zeroElement.desiredPosition.zero();
        m_data->desiredPositions.resize(newSize, zeroElement);
    }
};

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
    settingsStruct.staticTorquesCostActive = false;
    settingsStruct.forceMeanCostActive = true;
    settingsStruct.comCostActive = true;
    settingsStruct.comVelocityCostActive = true;
    settingsStruct.forceDerivativeCostActive = true;
    settingsStruct.pointAccelerationCostActive = true;
    settingsStruct.jointsRegularizationCostActive = true;
    settingsStruct.jointsVelocityCostActive = true;
    settingsStruct.swingCostActive = true;
    settingsStruct.phantomForcesCostActive = true;
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
    settingsStruct.comVelocityCostOverallWeight = 1.0;
    settingsStruct.comVelocityWeights(0) = 1.0;
    settingsStruct.comVelocityWeights(1) = 0.0;
    settingsStruct.comVelocityWeights(2) = 0.0;


    //    settingsStruct.minimumDt = 0.01;
    //    settingsStruct.controlPeriod = 0.1;
    //    settingsStruct.maximumDt = 0.05;

    settingsStruct.minimumDt = 0.1;
    settingsStruct.controlPeriod = 0.1;
    settingsStruct.maximumDt = 1.0;
    settingsStruct.horizon = 0.2;
    settingsStruct.activeControlPercentage = 1.0;

    settingsStruct.comCostActiveRange.setTimeInterval(settingsStruct.horizon*0, settingsStruct.horizon);

    iDynTree::VectorDynSize comPointReference(3);
    iDynTree::toEigen(comPointReference) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.1, 0.01, 0.0));
    auto comReference = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comPointReference);
    settingsStruct.desiredCoMTrajectory  = comReference;

    iDynTree::VectorDynSize comVelocityReference(3);
    iDynTree::toEigen(comVelocityReference) = iDynTree::toEigen(iDynTree::Position(0.1, 0.0, 0.0));
    auto comVelocityTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(comVelocityReference);
    settingsStruct.desiredCoMVelocityTrajectory  = comVelocityTrajectory;

    settingsStruct.meanPointPositionCostActiveRange.setTimeInterval(settingsStruct.horizon * 0, settingsStruct.horizon);
    MeanPointReferenceGenerator meanPointReferenceGenerator(2, 15.0);
    settingsStruct.desiredMeanPointPosition = meanPointReferenceGenerator.timeVaryingReference();
    settingsStruct.meanPointPositionCostTimeVaryingWeight = meanPointReferenceGenerator.timeVaryingWeight();
    iDynTree::toEigen(meanPointReferenceGenerator[0].desiredPosition) = iDynTree::toEigen(initialState.comPosition) + iDynTree::toEigen(iDynTree::Position(0.1, 0.0, 0.0));
    meanPointReferenceGenerator[0].desiredPosition(2) = 0.0;
    meanPointReferenceGenerator[0].activeRange.setTimeInterval(0.0, settingsStruct.horizon);
    meanPointReferenceGenerator[1].desiredPosition = meanPointReferenceGenerator[0].desiredPosition;
    meanPointReferenceGenerator[1].activeRange.setTimeInterval(settingsStruct.horizon + 1.0, settingsStruct.horizon + 1.0);



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
    iDynTree::toEigen(settingsStruct.velocityMaximumDerivative).setConstant(10);
    settingsStruct.velocityMaximumDerivative(0) = 10.0;
    settingsStruct.velocityMaximumDerivative(1) = 10.0;
    settingsStruct.planarVelocityHyperbolicTangentScaling = 0.1; //scales the position along z
    //    settingsStruct.normalVelocityHyperbolicSecantScaling = 1.0; //scales the force along z

    settingsStruct.complementarityDissipation = 10.0;

    settingsStruct.minimumCoMHeight = 0.5 * initialState.comPosition(2);

    ok = settings.setFromStruct(settingsStruct);
    ASSERT_IS_TRUE(ok);

    auto optimizerTest = std::make_shared<OptimizerTest>();


    ok = solver.setOptimizer(optimizerTest);
    ASSERT_IS_TRUE(ok);

//    auto integrator = std::make_shared<iDynTree::optimalcontrol::integrators::ForwardEuler>();

//    ok = solver.setIntegrator(integrator);
//    ASSERT_IS_TRUE(ok);

    ok = solver.specifySettings(settings);
    ASSERT_IS_TRUE(ok);

    ok = solver.setInitialState(initialState);
    ASSERT_IS_TRUE(ok);

    auto stateGuesses = std::make_shared<StateGuess>(comReference, initialState);
    auto controlGuesses = std::make_shared<DynamicalPlanner::TimeInvariantControl>(DynamicalPlanner::Control(vectorList.size(), settingsStruct.leftPointsPosition.size()));

    std::vector<DynamicalPlanner::State> optimalStates;
    std::vector<DynamicalPlanner::Control> optimalControls;

    for (size_t i = 0; i < 1; ++i) {
        ok = solver.solve(optimalStates, optimalControls);
        ASSERT_IS_TRUE(ok);
    }

    return EXIT_SUCCESS;
}
