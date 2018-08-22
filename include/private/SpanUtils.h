/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#ifndef DPLANNER_SPANUTILS_H
#define DPLANNER_SPANUTILS_H

#include <iDynTree/Core/Span.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <Eigen/Dense>

namespace DynamicalPlanner {
    namespace Private {

        inline Eigen::Map<Eigen::VectorXd> spanToEigen(iDynTree::Span<double> vec)
        {
            return Eigen::Map<Eigen::VectorXd>(vec.data(),vec.size());
        }

        inline Eigen::Map<const Eigen::VectorXd> spanToEigen(iDynTree::Span<const double> & vec)
        {
            return Eigen::Map<const Eigen::VectorXd>(vec.data(),vec.size());
        }

        typedef Eigen::Map<const Eigen::VectorXd> SpanEigenMap;

    }
}
#endif // DPLANNER_SPANUTILS_H
