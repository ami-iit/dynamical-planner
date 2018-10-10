/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/HyperbolicTangent.h>
#include <cmath>

using namespace DynamicalPlanner::Private;

HyperbolicTangent::HyperbolicTangent()
{ }

HyperbolicTangent::~HyperbolicTangent()
{ }

double HyperbolicTangent::eval(double x) const
{
    return std::tanh(m_K * x);
}

double HyperbolicTangent::evalDerivative(double x) const
{
    double sechKx = 1.0/std::cosh(m_K * x);
    return m_K * sechKx * sechKx;
}
