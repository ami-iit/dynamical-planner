/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_QUATERNIONUTILS_H
#define DPLANNER_QUATERNIONUTILS_H

#include <iDynTree/Core/MatrixFixSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/Rotation.h>

namespace DynamicalPlanner {
    namespace Private {

        bool QuaternionBoundsRespected(const iDynTree::Vector4 &quaternion);

        iDynTree::MatrixFixSize<4, 3> QuaternionLeftTrivializedDerivative(const iDynTree::Vector4& quaternion);

        iDynTree::MatrixFixSize<3, 4> QuaternionLeftTrivializedDerivativeInverse(const iDynTree::Vector4& quaternion);

        iDynTree::MatrixFixSize<4, 4> QuaternionLeftTrivializedDerivativeTimesOmegaJacobian(const iDynTree::Vector3& omega);

        iDynTree::MatrixFixSize<3, 4> QuaternionLeftTrivializedDerivativeInverseTimesQuaternionDerivativeJacobian(const iDynTree::Vector4& quatDerivative);

        iDynTree::Vector4 NormalizedQuaternion(const iDynTree::Vector4& quaternion);

        double QuaternionNorm(const iDynTree::Vector4 &quaternion);

        double QuaternionSquaredNorm(const iDynTree::Vector4& quaternion);

        iDynTree::Matrix4x4 NormalizedQuaternionDerivative(const iDynTree::Vector4& quaternion);

        iDynTree::MatrixFixSize<3,4> RotatedVectorQuaternionJacobian(const iDynTree::Vector3& originalVector, const iDynTree::Vector4& quaternion);

        iDynTree::Vector4 ErrorQuaternion(const iDynTree::Rotation& frameRotation, const iDynTree::Rotation& desiredRotation);

        iDynTree::Vector4 InverseQuaternion(const iDynTree::Vector4& quaternion);

        iDynTree::Matrix4x4 InverseQuaternionDerivative();

    }
}

#endif // DPLANNER_QUATERNIONUTILS_H
