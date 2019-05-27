/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SMOOTHINGFUNCTION_H
#define DPLANNER_SMOOTHINGFUNCTION_H

namespace DynamicalPlanner {
    namespace Private {
        class SmoothingFunction;
    }
}

class DynamicalPlanner::Private::SmoothingFunction {
protected:
    double m_K;
    bool m_disabled;
    double m_disabledValue;
public:
    SmoothingFunction();

    virtual ~SmoothingFunction();

    virtual void setScaling(double k);

    virtual double eval(double x) const = 0;

    virtual double evalDerivative(double x) const = 0;

    virtual double evalDoubleDerivative(double x) const = 0;

    virtual void disable(double constantValue = 1.0);
};


#endif // DPLANNER_SMOOTHINGFUNCTION_H
