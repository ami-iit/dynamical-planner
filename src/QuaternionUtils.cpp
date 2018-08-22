/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/QuaternionUtils.h>
#include <iDynTree/Core/EigenHelpers.h>

iDynTree::MatrixFixSize<4, 3> DynamicalPlanner::Private::QuaternionLeftTrivializedDerivative(iDynTree::Vector4 quaternion)
{
    iDynTree::MatrixFixSize<4, 3> outputMatrix;
    Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor> > map = iDynTree::toEigen(outputMatrix);
    map.topRows<1>() = -iDynTree::toEigen(quaternion).tail<3>().transpose();
    map.bottomRows<3>().setIdentity();
    map.bottomRows<3>() *= quaternion(0);
    map.bottomRows<3>() += iDynTree::skew(iDynTree::toEigen(quaternion).tail<3>());
    map *= 0.5;
    return outputMatrix;
}

iDynTree::MatrixFixSize<3, 4> DynamicalPlanner::Private::QuaternionLeftTrivializedDerivativeInverse(iDynTree::Vector4 quaternion)
{
    iDynTree::MatrixFixSize<3, 4> outputMatrix;
    Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> > map = iDynTree::toEigen(outputMatrix);

    map.setZero();
    map.leftCols<1>() = -iDynTree::toEigen(quaternion).tail<3>();
    map.rightCols<3>().setIdentity();
    map.rightCols<3>() *= iDynTree::toEigen(quaternion)(0);
    map.rightCols<3>() -= iDynTree::skew(iDynTree::toEigen(quaternion).tail<3>());

    map *= 2;
    return outputMatrix;
}

iDynTree::Vector4 DynamicalPlanner::Private::NormailizedQuaternion(iDynTree::Vector4 quaternion)
{
    iDynTree::Vector4 normalized;
    iDynTree::toEigen(normalized) = iDynTree::toEigen(quaternion).normalized();
    return normalized;
}

double DynamicalPlanner::Private::Norm(iDynTree::Vector4 quaternion)
{
    return iDynTree::toEigen(quaternion).norm();
}

double DynamicalPlanner::Private::SquaredNorm(iDynTree::Vector4 quaternion)
{
    return iDynTree::toEigen(quaternion).squaredNorm();
}

iDynTree::Matrix4x4 DynamicalPlanner::Private::NormalizedQuaternionDerivative(iDynTree::Vector4 quaternion)
{
    iDynTree::Matrix4x4 derivative;
    Eigen::Map<Eigen::Matrix<double, 4, 4, Eigen::RowMajor> > derivativeMap = iDynTree::toEigen(derivative);

    double threeHalfNorm = Norm(quaternion) * SquaredNorm(quaternion);
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> outerProduct = iDynTree::toEigen(quaternion) * iDynTree::toEigen(quaternion).transpose();
    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> identity;
    identity.setIdentity();
    derivativeMap = identity * SquaredNorm(quaternion) - outerProduct;
    derivativeMap *= 1/threeHalfNorm;
}
