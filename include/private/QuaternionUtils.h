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

namespace DynamicalPlanner {
    namespace Private {

        iDynTree::MatrixFixSize<4, 3> QuaternionLeftTrivializedDerivative(iDynTree::Vector4 quaternion);

        iDynTree::MatrixFixSize<3, 4> QuaternionLeftTrivializedDerivativeInverse(iDynTree::Vector4 quaternion);

        iDynTree::Vector4 NormailizedQuaternion(iDynTree::Vector4 quaternion);

        inline double Norm(iDynTree::Vector4 quaternion);

        inline double SquaredNorm(iDynTree::Vector4 quaternion);

        iDynTree::Matrix4x4 NormalizedQuaternionDerivative(iDynTree::Vector4 quaternion);
    }
}

#endif // DPLANNER_QUATERNIONUTILS_H
