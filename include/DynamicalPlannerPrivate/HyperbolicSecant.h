/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_HYPERBOLICSECANT_H
#define DPLANNER_HYPERBOLICSECANT_H

namespace DynamicalPlanner {
    namespace Private {
        class HyperbolicSecant;
    }
}

class DynamicalPlanner::Private::HyperbolicSecant {
    double m_K;
public:
    HyperbolicSecant();

    void setScaling(double k);

    double eval(double x) const;

    double evalDerivative(double x) const;
};

#endif // DPLANNER_HYPERBOLICSECANT_H
