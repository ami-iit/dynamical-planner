/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>

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
    iDynTree::getRandomVector(state.base_linearVelocity, -1.0, 1.0);
    iDynTree::getRandomVector(state.base_quaternionVelocity, -1.0, 1.0);
    iDynTree::getRandomVector(state.base_quaternion, -1.0, 1.0);
    state.base_quaternion(0) = std::abs(state.base_quaternion(0));
    return state;
}

void validateAdjoint(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    ExpressionsServer expressionServer(timelySharedKinDyn);
    expressionServer.updateRobotState(time, robotState);

    levi::Expression adjoint = expressionServer.adjointTransform("root_link", "l_sole");
    ASSERT_IS_TRUE(adjoint.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::Transform originalTransform = kinDyn->getRelativeTransform(robotState,"root_link", "l_sole");
    iDynTree::Matrix6x6 evalBuffer;
    iDynTree::toEigen(evalBuffer) = adjoint.evaluate();
    ASSERT_EQUAL_MATRIX(evalBuffer, originalTransform.asAdjointTransform());

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation;
    iDynTree::Vector6 perturbedCol, firstOrderTaylor;
    Eigen::MatrixXd derivative;

    originalJoints = iDynTree::toEigen(robotState.s);

    for (Eigen::Index col = 0; col < 6; ++col) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        expressionServer.updateRobotState(time, robotState);

        derivative = adjoint.getColumnDerivative(col, expressionServer.jointsPosition()).evaluate();

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
    ExpressionsServer expressionServer(timelySharedKinDyn);
    expressionServer.updateRobotState(time, robotState);

    levi::Expression adjoint = expressionServer.adjointTransformWrench("root_link", "l_sole");
    ASSERT_IS_TRUE(adjoint.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::Transform originalTransform = kinDyn->getRelativeTransform(robotState,"root_link", "l_sole");
    iDynTree::Matrix6x6 evalBuffer;
    iDynTree::toEigen(evalBuffer) = adjoint.evaluate();
    ASSERT_EQUAL_MATRIX(evalBuffer, originalTransform.asAdjointTransformWrench());

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation;
    iDynTree::Vector6 perturbedCol, firstOrderTaylor;
    Eigen::MatrixXd derivative;

    originalJoints = iDynTree::toEigen(robotState.s);

    for (Eigen::Index col = 0; col < 6; ++col) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        expressionServer.updateRobotState(time, robotState);

        derivative = adjoint.getColumnDerivative(col, expressionServer.jointsPosition()).evaluate();

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
    ExpressionsServer expressionServer(timelySharedKinDyn);
    expressionServer.updateRobotState(time, robotState);

    levi::Expression jacobian = expressionServer.relativeLeftJacobian("root_link", "l_sole");
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
        expressionServer.updateRobotState(time, robotState);


        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        derivative = jacobian.getColumnDerivative(col, expressionServer.jointsPosition()).evaluate();
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
    ExpressionsServer expressionServer(timelySharedKinDyn);
    expressionServer.updateRobotState(time, robotState);

    TransformExpression newTransform = expressionServer.relativeTransform("root_link", "l_sole");

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

    jacobian = transformedPos.getColumnDerivative(0, expressionServer.jointsPosition()).evaluate();

    originalJoints = iDynTree::toEigen(robotState.s);

    for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        robotState.s(joint) += perturbation;

        perturbedPosition = kinDyn->getRelativeTransform(robotState, "root_link", "l_sole") * testPosition;

        firstOrderTaylor = iDynTree::toEigen(testTransformedPos) + jacobian * (iDynTree::toEigen(robotState.s) - originalJoints);

        ASSERT_EQUAL_VECTOR_TOL(perturbedPosition, firstOrderTaylor, perturbation/1000.0);

    }
}

void validateCom(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    ExpressionsServer expressionServer(timelySharedKinDyn);
    expressionServer.updateRobotState(time, robotState);

    levi::Expression comPosition = expressionServer.worldToBase() * expressionServer.comInBase();

    ASSERT_IS_TRUE(comPosition.isValidExpression());

    levi::Expression comJacobian = comPosition.getColumnDerivative(0, expressionServer.jointsPosition());

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
        jacobianDerivative = comJacobian.getColumnDerivative(col, expressionServer.jointsPosition()).evaluate();
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

void validateRelativeVelocityExpression(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());
    ExpressionsServer expressionServer(timelySharedKinDyn);
    expressionServer.updateRobotState(time, robotState);

    levi::Expression relativeVelocityExpression = expressionServer.relativeVelocity("root_link", "r_sole");

    ASSERT_IS_TRUE(relativeVelocityExpression.isValidExpression());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    iDynTree::LinVelocity baseLinVelocity, baseAngVelocity;
    baseLinVelocity = robotState.base_linearVelocity;
    iDynTree::Vector4 quaternionNormalized;
    iDynTree::toEigen(quaternionNormalized) = iDynTree::toEigen(robotState.base_quaternion).normalized();
    iDynTree::toEigen(baseAngVelocity) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(quaternionNormalized)) *
        iDynTree::toEigen(robotState.base_quaternionVelocity);

    iDynTree::Twist checkRelativeVelocity = kinDyn->getFrameVel(robotState, "r_sole", iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION) -
        kinDyn->getRelativeTransform(robotState, "r_sole", "root_link") * iDynTree::Twist(baseLinVelocity, baseAngVelocity);

    iDynTree::Vector6 relVelocityEval;
    iDynTree::toEigen(relVelocityEval) = relativeVelocityExpression.evaluate();
    ASSERT_EQUAL_VECTOR(relVelocityEval, checkRelativeVelocity);

    iDynTree::toEigen(relVelocityEval) = relativeVelocityExpression.getColumnDerivative(0, expressionServer.jointsVelocity()).evaluate() * expressionServer.jointsVelocity().evaluate();
    ASSERT_EQUAL_VECTOR(relVelocityEval, checkRelativeVelocity);

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation, firstOrderTaylor(3);
    iDynTree::Twist perturbedVelocity;

    Eigen::MatrixXd jacobian(6, robotState.s.size());

    jacobian = relativeVelocityExpression.getColumnDerivative(0, expressionServer.jointsPosition()).evaluate();

    originalJoints = iDynTree::toEigen(robotState.s);

    for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        robotState.s(joint) += perturbation;

        perturbedVelocity = kinDyn->getFrameVel(robotState, "r_sole", iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION) -
            kinDyn->getRelativeTransform(robotState, "r_sole", "root_link") * iDynTree::Twist(baseLinVelocity, baseAngVelocity);

        firstOrderTaylor = iDynTree::toEigen(checkRelativeVelocity) + jacobian * (iDynTree::toEigen(robotState.s) - originalJoints);

        ASSERT_EQUAL_VECTOR_TOL(perturbedVelocity, firstOrderTaylor, perturbation/1000.0);

    }
}

void validateQuaternionError(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, double time) {
    std::shared_ptr<ExpressionsServer> expressionsServer = std::make_shared<ExpressionsServer>(timelySharedKinDyn);
    RobotState robotState = RandomRobotState(timelySharedKinDyn->model());

    expressionsServer->updateRobotState(time, robotState);

    iDynTree::Rotation desiredRotation = iDynTree::getRandomRotation();
    levi::Variable quat_des(4, "quat_des");
    quat_des = iDynTree::toEigen(desiredRotation.asQuaternion());

    SharedKinDynComputationsPointer kinDyn = timelySharedKinDyn->get(time);

    levi::Expression quaternionErrorExpression = expressionsServer->quaternionError("neck_2", quat_des);

    ASSERT_IS_TRUE(quaternionErrorExpression.isValidExpression());

    iDynTree::Transform frameTransform = kinDyn->getWorldTransform(robotState, "neck_2");

    iDynTree::Vector4 quaternion_iDyn = ErrorQuaternion(frameTransform.getRotation(), desiredRotation);

    iDynTree::Vector4 quaternion_levi;

    iDynTree::toEigen(quaternion_levi) = quaternionErrorExpression.evaluate();

    ASSERT_EQUAL_VECTOR(quaternion_levi, quaternion_iDyn);

    double perturbation = 1e-3;
    Eigen::VectorXd originalJoints, jointsPerturbation;
    iDynTree::Vector4 firstOrderTaylor, perturbedQuaternion, originalQuaternion;
    Eigen::MatrixXd derivative;

    derivative = quaternionErrorExpression.getColumnDerivative(0, expressionsServer->jointsPosition()).evaluate();

    originalJoints = iDynTree::toEigen(robotState.s);

    for (unsigned int joint = 0; joint < robotState.s.size(); ++joint) {

        iDynTree::toEigen(robotState.s) = originalJoints;
        robotState.s(joint) += perturbation;

        frameTransform = kinDyn->getWorldTransform(robotState, "neck_2");
        perturbedQuaternion = ErrorQuaternion(frameTransform.getRotation(), desiredRotation);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(quaternion_iDyn) + derivative * (iDynTree::toEigen(robotState.s) - originalJoints);

        ASSERT_EQUAL_VECTOR_TOL(perturbedQuaternion, firstOrderTaylor, perturbation/1000.0);

    }

    iDynTree::toEigen(robotState.s) = originalJoints;

    Eigen::MatrixXd derivativeBase = quaternionErrorExpression.getColumnDerivative(0, expressionsServer->baseQuaternion()).evaluate();

    originalQuaternion = robotState.base_quaternion;

    for (unsigned int i = 0; i < 4; ++i) {

        robotState.base_quaternion = originalQuaternion;
        robotState.base_quaternion(i) += perturbation;

        frameTransform = kinDyn->getWorldTransform(robotState, "neck_2");
        perturbedQuaternion = ErrorQuaternion(frameTransform.getRotation(), desiredRotation);

        iDynTree::toEigen(firstOrderTaylor) = iDynTree::toEigen(quaternion_iDyn) + derivativeBase * (iDynTree::toEigen(robotState.base_quaternion) - iDynTree::toEigen(originalQuaternion));

        ASSERT_EQUAL_VECTOR_TOL(perturbedQuaternion, firstOrderTaylor, perturbation/1000.0);

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

    validateRelativeVelocityExpression(timelySharedKinDyn, 0.0);

    validateRelativeVelocityExpression(timelySharedKinDyn, 1.0);

    validateQuaternionError(timelySharedKinDyn, 0.0);

    return 0;
}
