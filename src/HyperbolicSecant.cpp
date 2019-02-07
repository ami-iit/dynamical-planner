/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/HyperbolicSecant.h>
#include <cmath>

using namespace DynamicalPlanner::Private;

HyperbolicSecant::HyperbolicSecant()
{ }

HyperbolicSecant::~HyperbolicSecant()
{ }

double HyperbolicSecant::eval(double x) const
{
    return 1.0/std::cosh(m_K * x);
}

double HyperbolicSecant::evalDerivative(double x) const
{
    return -eval(x) * std::tanh(m_K * x) * m_K;
}

double HyperbolicSecant::evalDoubleDerivative(double x) const
{
    double tanhKx = std::tanh(m_K * x);
    return -evalDerivative(x) * tanhKx * m_K - eval(x) * (1 - tanhKx * tanhKx) * m_K * m_K;
}
