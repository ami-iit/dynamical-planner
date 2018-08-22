/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/QuaternionUtils.h>

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

int main() {
    validateQuaternionLeftTrivializedDerivative(iDynTree::getRandomRotation());
    return EXIT_SUCCESS;
}

