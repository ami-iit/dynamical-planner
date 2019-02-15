/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeJacobianExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/TransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/CoMInBaseExpression.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>

#include <URDFdir.h>

#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/Rotation.h>
#include <iDynTree/ModelIO/ModelLoader.h>

#include <cmath>
#include <chrono>
#include <iostream>

using namespace DynamicalPlanner::Private;


void validateQuaternionExpressions(const iDynTree::Rotation & R) {

    levi::Variable q(4, "q");

    iDynTree::Vector4 quaternion = R.asQuaternion();
    q = iDynTree::toEigen(quaternion);
    ASSERT_EQUAL_MATRIX(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(quaternion),
                        (2.0 * (DynamicalPlanner::Private::E_Expression(q))).evaluate());

    ASSERT_EQUAL_MATRIX(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverse(quaternion),
                        (2.0 * (DynamicalPlanner::Private::G_Expression(q))).evaluate());
}

void configureSharedKinDyn(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn) {
    std::vector<std::string> vectorList({"torso_pitch", "torso_roll", "torso_yaw", "l_shoulder_pitch", "l_shoulder_roll",
                                         "l_shoulder_yaw", "l_elbow", "r_shoulder_pitch", "r_shoulder_roll", "r_shoulder_yaw",
                                         "r_elbow", "l_hip_pitch", "l_hip_roll", "l_hip_yaw", "l_knee", "l_ankle_pitch",
                                         "l_ankle_roll", "r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

    //    std::vector<std::string> vectorList({"r_hip_pitch", "r_hip_roll", "r_hip_yaw", "r_knee", "r_ankle_pitch", "r_ankle_roll"});

    iDynTree::ModelLoader modelLoader;
    bool ok = modelLoader.loadModelFromFile(getAbsModelPath("iCubGenova04.urdf"));
    ASSERT_IS_TRUE(ok);
    ok = modelLoader.loadReducedModelFromFullModel(modelLoader.model(), vectorList);
    ASSERT_IS_TRUE(ok);
    assert(timelySharedKinDyn);
    ok = timelySharedKinDyn->loadRobotModel(modelLoader.model());
    ASSERT_IS_TRUE(ok);
    //    ASSERT_IS_TRUE(sharedKinDyn->model().getNrOfDOFs() == 23);
    std::vector<double> timings(2);
    timings[0] = 0.0;
    timings[1] = 1.0;

    ok = timelySharedKinDyn->setTimings(timings);
    ASSERT_IS_TRUE(ok);
}

RobotState RandomRobotState(const iDynTree::Model& model) {
    RobotState state;
    state.s.resize(static_cast<unsigned int>(model.getNrOfJoints()));
    iDynTree::getRandomVector(state.s, -1.0, 1.0);
    state.s_dot.resize(static_cast<unsigned int>(model.getNrOfJoints()));
    iDynTree::getRandomVector(state.s_dot, -1.0, 1.0);
    state.base_position = iDynTree::getRandomPosition();
    state.base_velocity = iDynTree::getRandomTwist();
    iDynTree::getRandomVector(state.base_quaternion, -1.0, 1.0);
    state.base_quaternion(0) = std::abs(state.base_quaternion(0));
    return state;
}

void validateAdjoint(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    levi::Variable q(timelySharedKinDyn->model().getNrOfJoints(), "q");
    levi::ScalarVariable t("t");

    t = time;
    q = iDynTree::toEigen(robotState.s);

    levi::Expression adjoint = AdjointTransformExpression(timelySharedKinDyn, &robotState, "root_link", "l_sole", q, t);
    ASSERT_IS_TRUE(adjoint.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::Transform originalTransform = kinDyn->getRelativeTransform(robotState,"root_link", "l_sole");
    ASSERT_IS_TRUE(adjoint.evaluate() == iDynTree::toEigen(originalTransform.asAdjointTransform()));

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation;
    iDynTree::Vector6 perturbedCol, firstOrderTaylor;
    Eigen::MatrixXd derivative;

    originalJoints = iDynTree::toEigen(robotState.s);

    for (Eigen::Index col = 0; col < 6; ++col) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        q = originalJoints;

        derivative = adjoint.getColumnDerivative(col, q).evaluate();

        for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

            iDynTree::toEigen(robotState.s) = originalJoints;
            robotState.s(joint) += perturbation;

            iDynTree::toEigen(perturbedCol) = iDynTree::toEigen(kinDyn->getRelativeTransform(robotState,"root_link", "l_sole").asAdjointTransform()).col(col);

            iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalTransform.asAdjointTransform()).col(col) + derivative * (iDynTree::toEigen(robotState.s) - originalJoints);

            ASSERT_EQUAL_VECTOR_TOL(perturbedCol, firstOrderTaylor, perturbation/10.0);

        }
    }
}

void validateAdjointWrench(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    levi::Variable q(timelySharedKinDyn->model().getNrOfJoints(), "q");
    levi::ScalarVariable t("t");

    t = time;
    q = iDynTree::toEigen(robotState.s);

    levi::Expression adjoint = AdjointTransformWrenchExpression(timelySharedKinDyn, &robotState, "root_link", "l_sole", q, t);
    ASSERT_IS_TRUE(adjoint.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::Transform originalTransform = kinDyn->getRelativeTransform(robotState,"root_link", "l_sole");
    ASSERT_IS_TRUE(adjoint.evaluate() == iDynTree::toEigen(originalTransform.asAdjointTransformWrench()));

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation;
    iDynTree::Vector6 perturbedCol, firstOrderTaylor;
    Eigen::MatrixXd derivative;

    originalJoints = iDynTree::toEigen(robotState.s);

    for (Eigen::Index col = 0; col < 6; ++col) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        q = originalJoints;

        derivative = adjoint.getColumnDerivative(col, q).evaluate();

        for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

            iDynTree::toEigen(robotState.s) = originalJoints;
            robotState.s(joint) += perturbation;

            iDynTree::toEigen(perturbedCol) = iDynTree::toEigen(kinDyn->getRelativeTransform(robotState,"root_link", "l_sole").asAdjointTransformWrench()).col(col);

            iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalTransform.asAdjointTransformWrench()).col(col) + derivative * (iDynTree::toEigen(robotState.s) - originalJoints);

            ASSERT_EQUAL_VECTOR_TOL(perturbedCol, firstOrderTaylor, perturbation/10.0);

        }
    }
}

void validateJacobian(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    levi::Variable q(timelySharedKinDyn->model().getNrOfJoints(), "q");
    levi::ScalarVariable t("t");

    t = time;
    q = iDynTree::toEigen(robotState.s);

    levi::Expression jacobian = RelativeLeftJacobianExpression(timelySharedKinDyn, &robotState, "root_link", "l_sole", q, t);
    ASSERT_IS_TRUE(jacobian.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::MatrixDynSize originalJacobian(6, robotState.s.size()), perturbedJacobian = originalJacobian;
    iDynTree::FrameIndex baseFrame = timelySharedKinDyn->model().getFrameIndex("root_link");
    iDynTree::FrameIndex targetFrame = timelySharedKinDyn->model().getFrameIndex("l_sole");

    bool ok = kinDyn->getRelativeJacobian(robotState, baseFrame, targetFrame, originalJacobian, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
    ASSERT_IS_TRUE(ok);
    ASSERT_IS_TRUE(jacobian.evaluate() == iDynTree::toEigen(originalJacobian));

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation;
    iDynTree::Vector6 perturbedCol, firstOrderTaylor;
    Eigen::MatrixXd derivative;

    originalJoints = iDynTree::toEigen(robotState.s);

    for (Eigen::Index col = 0; col < robotState.s.size(); ++col) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        q = originalJoints;


        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        derivative = jacobian.getColumnDerivative(col, q).evaluate();
        std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (col " << col << "): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

        for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

            iDynTree::toEigen(robotState.s) = originalJoints;
            robotState.s(joint) += perturbation;

            ok = kinDyn->getRelativeJacobian(robotState, baseFrame, targetFrame, perturbedJacobian, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
            ASSERT_IS_TRUE(ok);
            iDynTree::toEigen(perturbedCol) = iDynTree::toEigen(perturbedJacobian).col(col);

            iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(originalJacobian).col(col) + derivative * (iDynTree::toEigen(robotState.s) - originalJoints);

            ASSERT_EQUAL_VECTOR_TOL(perturbedCol, firstOrderTaylor, perturbation/100.0);

        }
    }

}

void validateTransform(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    levi::Variable q(timelySharedKinDyn->model().getNrOfJoints(), "q");
    levi::ScalarVariable t("t");

    t = time;
    q = iDynTree::toEigen(robotState.s);

    TransformExpression newTransform = TransformExpression::RelativeTransform(timelySharedKinDyn, &robotState, "root_link", "l_sole", q, t);

    iDynTree::Position testPosition = iDynTree::getRandomPosition();

    levi::Constant testExpr(3,1,"pos");
    testExpr = iDynTree::toEigen(testPosition);

    levi::Expression transformedPos = newTransform * testExpr;

    ASSERT_IS_TRUE(transformedPos.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::Position testTransformedPos = kinDyn->getRelativeTransform(robotState, "root_link", "l_sole") * testPosition;

    iDynTree::Position positionFromExpression;

    iDynTree::toEigen(positionFromExpression) = transformedPos.evaluate();

    ASSERT_EQUAL_VECTOR(testTransformedPos, positionFromExpression);

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation, firstOrderTaylor(3);
    iDynTree::Position perturbedPosition;

    Eigen::MatrixXd jacobian(3, robotState.s.size());

    jacobian = transformedPos.getColumnDerivative(0, q).evaluate();

    originalJoints = iDynTree::toEigen(robotState.s);

    for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        robotState.s(joint) += perturbation;

        perturbedPosition = kinDyn->getRelativeTransform(robotState, "root_link", "l_sole") * testPosition;

        firstOrderTaylor = iDynTree::toEigen(testTransformedPos) + jacobian * (iDynTree::toEigen(robotState.s) - originalJoints);

        ASSERT_EQUAL_VECTOR_TOL(perturbedPosition, firstOrderTaylor, perturbation/10.0);

    }


}

void validateCom(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    levi::Variable q(timelySharedKinDyn->model().getNrOfJoints(), "q");
    levi::ScalarVariable t("t");

    t = time;
    q = iDynTree::toEigen(robotState.s);

    levi::Variable p_base(3, "p_b");

    p_base = iDynTree::toEigen(robotState.base_position);

    levi::Variable base_quat(4, "qb");

    base_quat = iDynTree::toEigen(robotState.base_quaternion);

    levi::Expression baseRotation = RotationExpression(NormalizedQuaternion(base_quat).asVariable());

    TransformExpression baseTransform(p_base, baseRotation);

    levi::Expression comPosition = baseTransform * CoMInBaseExpression(timelySharedKinDyn, &robotState, q, t);

    ASSERT_IS_TRUE(comPosition.isValidExpression());

    levi::Expression comJacobian = comPosition.getColumnDerivative(0, q);

    iDynTree::MatrixDynSize iDyntree_jacobian;

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    bool ok = kinDyn->getCenterOfMassJacobian(robotState, iDyntree_jacobian);
    ASSERT_IS_TRUE(ok);

    iDynTree::MatrixDynSize jointsJacobian(3, robotState.s.size());

    iDynTree::toEigen(jointsJacobian) = iDynTree::toEigen(iDyntree_jacobian).block(0, 6, 3, robotState.s.size());

    iDynTree::MatrixDynSize exprJacobian = jointsJacobian;
    exprJacobian.zero();

    iDynTree::toEigen(exprJacobian) = comJacobian.evaluate();

    ASSERT_EQUAL_MATRIX(jointsJacobian, exprJacobian);

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation, firstOrderTaylor(3);
    iDynTree::VectorDynSize perturbedCol(3);
    Eigen::MatrixXd originalJacobian = iDynTree::toEigen(iDyntree_jacobian);

    originalJoints = iDynTree::toEigen(robotState.s);

    for (unsigned int col = 0; col < robotState.s.size(); ++col) {
        Eigen::MatrixXd jacobianDerivative(3, robotState.s.size());

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        jacobianDerivative = comJacobian.getColumnDerivative(col, q).evaluate();
        std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (Com jacobian col " << col << "): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

        for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

            iDynTree::toEigen(robotState.s) = originalJoints;
            robotState.s(joint) += perturbation;

            ok = kinDyn->getCenterOfMassJacobian(robotState, iDyntree_jacobian);
            ASSERT_IS_TRUE(ok);

            iDynTree::toEigen(perturbedCol) = iDynTree::toEigen(iDyntree_jacobian).col(col+6);

            firstOrderTaylor = originalJacobian.col(col+6) + jacobianDerivative * (iDynTree::toEigen(robotState.s) - originalJoints);

            ASSERT_EQUAL_VECTOR_TOL(perturbedCol, firstOrderTaylor, perturbation/10.0);

        }
    }


}


int main() {

    std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn = std::make_shared<TimelySharedKinDynComputations>();
    configureSharedKinDyn(timelySharedKinDyn);

    validateQuaternionExpressions(iDynTree::getRandomRotation());

    validateAdjoint(timelySharedKinDyn, 0.0);

    validateAdjoint(timelySharedKinDyn, 1.0);

    validateAdjointWrench(timelySharedKinDyn, 0.0);

    validateAdjointWrench(timelySharedKinDyn, 1.0);

    validateJacobian(timelySharedKinDyn, 0.0);

    validateJacobian(timelySharedKinDyn, 1.0);

    validateTransform(timelySharedKinDyn, 0.0);

    validateTransform(timelySharedKinDyn, 1.0);

    validateCom(timelySharedKinDyn, 0.0);

    validateCom(timelySharedKinDyn, 1.0);

    return 0;
}
