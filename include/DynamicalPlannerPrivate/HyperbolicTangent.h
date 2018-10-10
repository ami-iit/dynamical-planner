/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_HYPERBOLICTANGENT_H
#define DPLANNER_HYPERBOLICTANGENT_H

#include <DynamicalPlannerPrivate/SmoothingFunction.h>

namespace DynamicalPlanner {
    namespace Private {
        class HyperbolicTangent;
    }
}

class DynamicalPlanner::Private::HyperbolicTangent : public DynamicalPlanner::Private::SmoothingFunction {
public:
    HyperbolicTangent();

    virtual ~HyperbolicTangent() final;

    virtual double eval(double x) const final;

    virtual double evalDerivative(double x) const final;
};

#endif // DPLANNER_HYPERBOLICTANGENT_H
