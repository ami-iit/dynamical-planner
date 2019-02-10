/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_HYPERBOLICSECANT_H
#define DPLANNER_HYPERBOLICSECANT_H

#include <DynamicalPlannerPrivate/Utilities/SmoothingFunction.h>

namespace DynamicalPlanner {
    namespace Private {
        class HyperbolicSecant;
    }
}

class DynamicalPlanner::Private::HyperbolicSecant : public DynamicalPlanner::Private::SmoothingFunction {
public:
    HyperbolicSecant();

    virtual ~HyperbolicSecant() final;

    virtual double eval(double x) const final;

    virtual double evalDerivative(double x) const final;

    virtual double evalDoubleDerivative(double x) const final;
};

#endif // DPLANNER_HYPERBOLICSECANT_H
