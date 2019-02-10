/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>

#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/Rotation.h>

void validateQuaternionLeftTrivializedDerivative(const iDynTree::Rotation & R)
{
    double numericalDerivStep = 1e-8;
    iDynTree::Vector4 quaternion = R.asQuaternion();
    iDynTree::Vector4 quaternionPerturbedUp;
    iDynTree::Vector4 quaternionPerturbedDown;

    iDynTree::MatrixFixSize<4, 3> OmegaToDotQuat = DynamicalPlanner::Private::QuaternionLeftTrivializedDerivative(quaternion);
    iDynTree::MatrixFixSize<3, 4> DotQuatToOmega = DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverse(quaternion);

    iDynTree::Matrix3x3 identityCheck;
    iDynTree::toEigen(identityCheck) = iDynTree::toEigen(DotQuatToOmega) * iDynTree::toEigen(OmegaToDotQuat);

    ASSERT_EQUAL_MATRIX(identityCheck, iDynTree::Rotation::Identity());

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbedUp = quaternion;
        quaternionPerturbedUp(i) = quaternionPerturbedUp(i) + numericalDerivStep;
        quaternionPerturbedDown = quaternion;
        quaternionPerturbedDown(i) = quaternionPerturbedDown(i) - numericalDerivStep;
        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
        toEigen(quaternionPerturbedUp) = toEigen(quaternionPerturbedUp).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbedUp) = toEigen(quaternionPerturbedUp).array().max(minQuaternion.array());
        toEigen(quaternionPerturbedUp).normalize();

        toEigen(quaternionPerturbedDown) = toEigen(quaternionPerturbedDown).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbedDown) = toEigen(quaternionPerturbedDown).array().max(minQuaternion.array());
        toEigen(quaternionPerturbedDown).normalize();


        iDynTree::Vector4 quatNumDev;
        toEigen(quatNumDev) = (toEigen(quaternionPerturbedUp)-toEigen(quaternionPerturbedDown))/(2*numericalDerivStep);

        iDynTree::Rotation RPerturbedUp = iDynTree::Rotation::RotationFromQuaternion(quaternionPerturbedUp);
        iDynTree::Rotation RPerturbedDown = iDynTree::Rotation::RotationFromQuaternion(quaternionPerturbedDown);

        iDynTree::Matrix3x3 RNumDev;
        toEigen(RNumDev) = (toEigen(RPerturbedUp)-toEigen(RPerturbedDown))/(2*numericalDerivStep);

        iDynTree::Matrix3x3 skewMat;
        toEigen(skewMat) = toEigen(R).transpose()*toEigen(RNumDev);

        iDynTree::Vector3 angVelNum;
        toEigen(angVelNum) = iDynTree::unskew(toEigen(skewMat));

        iDynTree::Vector4 checkDotQuat;
        toEigen(checkDotQuat) = toEigen(OmegaToDotQuat) * toEigen(angVelNum);

        ASSERT_EQUAL_VECTOR_TOL(checkDotQuat,quatNumDev,numericalDerivStep * 10);

        iDynTree::Vector3 angVelCheck;
        toEigen(angVelCheck) = toEigen(DotQuatToOmega) * toEigen(quatNumDev);

        ASSERT_EQUAL_VECTOR_TOL(angVelCheck,angVelNum,numericalDerivStep * 10);

    }
}

void validateMapJacobian(const iDynTree::Rotation & R) {
    double perturbationValue = 1e-2;
    iDynTree::Vector4 quaternion = R.asQuaternion();
    iDynTree::Vector4 quaternionPerturbed, perturbation, firstOrderTaylor, output, perturbedOutput;
    iDynTree::Vector3 omega;
    iDynTree::getRandomVector(omega);

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    toEigen(output) = toEigen(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivative(quaternion)) * toEigen(omega);


    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbed = quaternion;
        quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().max(minQuaternion.array());
        toEigen(quaternionPerturbed).normalize();

        toEigen(perturbation) = toEigen(quaternionPerturbed) - toEigen(quaternion);

        toEigen(perturbedOutput) = toEigen(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivative(quaternionPerturbed)) * toEigen(omega);

        toEigen(firstOrderTaylor) = toEigen(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeTimesOmegaJacobian(omega)) * toEigen(perturbation);

        toEigen(firstOrderTaylor) += toEigen(output);

        ASSERT_EQUAL_VECTOR_TOL(perturbedOutput,firstOrderTaylor,perturbationValue/10);

    }

}

void validateInverseMapJacobian(const iDynTree::Rotation & R) {
    double perturbationValue = 1e-2;
    iDynTree::Vector4 quaternion = R.asQuaternion();
    iDynTree::Vector4 quaternionPerturbed, perturbation, quaternionDerivative;
    iDynTree::Vector3 omega, omegaPerturbed, firstOrderTaylor;
    iDynTree::getRandomVector(quaternionDerivative);

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    toEigen(omega) = toEigen(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverse(quaternion)) * toEigen(quaternionDerivative);


    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbed = quaternion;
        quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().max(minQuaternion.array());
        toEigen(quaternionPerturbed).normalize();

        toEigen(perturbation) = toEigen(quaternionPerturbed) - toEigen(quaternion);

        toEigen(omegaPerturbed) = toEigen(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverse(quaternionPerturbed)) *
                toEigen(quaternionDerivative);

        toEigen(firstOrderTaylor) =
                toEigen(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverseTimesQuaternionDerivativeJacobian(quaternionDerivative)) *
                toEigen(perturbation);

        toEigen(firstOrderTaylor) += toEigen(omega);

        ASSERT_EQUAL_VECTOR_TOL(omegaPerturbed,firstOrderTaylor,perturbationValue/100);

    }
}

void validateNormalizeQuaternionJacobian() {
    double perturbationValue = 1e-2;
    iDynTree::Vector4 quaternion;
    iDynTree::getRandomVector(quaternion, -1.0, 1.0);
    quaternion(0) = iDynTree::getRandomDouble();

    iDynTree::Vector4 originalQuaternionNormalized = DynamicalPlanner::Private::NormalizedQuaternion(quaternion);

    iDynTree::Vector4 quaternionPerturbed, perturbation, firstOrderTaylor, perturbedQuaternionNormalized;

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbed = quaternion;
        quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().max(minQuaternion.array());

        toEigen(perturbation) = toEigen(quaternionPerturbed) - toEigen(quaternion);

        perturbedQuaternionNormalized = DynamicalPlanner::Private::NormalizedQuaternion(quaternionPerturbed);

        toEigen(firstOrderTaylor) = toEigen(DynamicalPlanner::Private::NormalizedQuaternionDerivative(quaternion)) * toEigen(perturbation);

        toEigen(firstOrderTaylor) += toEigen(originalQuaternionNormalized);

        ASSERT_EQUAL_VECTOR_TOL(perturbedQuaternionNormalized, firstOrderTaylor,perturbationValue/10);

    }

}

void validateRotatedVectorJacobian(const iDynTree::Rotation & R) {
    double perturbationValue = 1e-3;
    iDynTree::Vector4 quaternion = R.asQuaternion();
    iDynTree::Vector4 quaternionPerturbed, perturbation;
    iDynTree::Rotation perturbedRotation;
    iDynTree::Vector3 vector, output, perturbedOutput, firstOrderTaylor;
    iDynTree::getRandomVector(vector, -10, 10);

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    toEigen(output) = toEigen(R) * toEigen(vector);


    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbed = quaternion;
        quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().max(minQuaternion.array());
        toEigen(quaternionPerturbed).normalize();

        toEigen(perturbation) = toEigen(quaternionPerturbed) - toEigen(quaternion);

        perturbedRotation.fromQuaternion(quaternionPerturbed);

        toEigen(perturbedOutput) = toEigen(perturbedRotation) * toEigen(vector);

        toEigen(firstOrderTaylor) = toEigen(DynamicalPlanner::Private::RotatedVectorQuaternionJacobian(vector, quaternion)) * toEigen(perturbation);

        toEigen(firstOrderTaylor) += toEigen(output);

        ASSERT_EQUAL_VECTOR_TOL(perturbedOutput,firstOrderTaylor,perturbationValue/10);

    }

}

void validateInverseQuaternion(const iDynTree::Rotation & R) {
    double perturbationValue = 1e-3;
    iDynTree::Vector4 quaternion = R.asQuaternion();
    iDynTree::Vector4 quaternionPerturbed, perturbation;
    iDynTree::Vector4 invertedQuaternion, identityCheck, perturbedInvertedQuaternion, firstOrderTaylor;

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    invertedQuaternion = DynamicalPlanner::Private::InverseQuaternion(quaternion);


    identityCheck = (iDynTree::Rotation::RotationFromQuaternion(quaternion) *
                     iDynTree::Rotation::RotationFromQuaternion(invertedQuaternion).inverse()).asQuaternion();

    ASSERT_EQUAL_VECTOR(identityCheck, iDynTree::Rotation::Identity().asQuaternion());



    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbed = quaternion;
        quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().min(maxQuaternion.array());
        toEigen(quaternionPerturbed) = toEigen(quaternionPerturbed).array().max(minQuaternion.array());
        toEigen(quaternionPerturbed).normalize();

        toEigen(perturbation) = toEigen(quaternionPerturbed) - toEigen(quaternion);

        perturbedInvertedQuaternion = DynamicalPlanner::Private::InverseQuaternion(quaternionPerturbed);

        toEigen(firstOrderTaylor) = toEigen(DynamicalPlanner::Private::InverseQuaternionDerivative()) * toEigen(perturbation);

        toEigen(firstOrderTaylor) += toEigen(invertedQuaternion);

        ASSERT_EQUAL_VECTOR_TOL(perturbedInvertedQuaternion,firstOrderTaylor,perturbationValue/10);

    }

}

void validateExpressions(const iDynTree::Rotation & R) {

    levi::Variable q(4, "q");
    DynamicalPlanner::Private::E_Expression E_levi(q);
    DynamicalPlanner::Private::G_Expression G_levi(q);


    iDynTree::Vector4 quaternion = R.asQuaternion();
    q = iDynTree::toEigen(quaternion);
    ASSERT_EQUAL_MATRIX(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(quaternion),
                        (2.0 * (*E_levi)).evaluate());

    ASSERT_EQUAL_MATRIX(DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverse(quaternion),
                        (2.0 * (*G_levi)).evaluate());
}

int main() {
    validateQuaternionLeftTrivializedDerivative(iDynTree::getRandomRotation());
    validateMapJacobian(iDynTree::getRandomRotation());
    validateInverseMapJacobian(iDynTree::getRandomRotation());
    validateNormalizeQuaternionJacobian();
    validateRotatedVectorJacobian(iDynTree::getRandomRotation());
    validateExpressions(iDynTree::getRandomRotation());

    return EXIT_SUCCESS;
}

