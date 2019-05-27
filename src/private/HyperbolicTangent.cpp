/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Utilities/HyperbolicTangent.h>
#include <cmath>

using namespace DynamicalPlanner::Private;

HyperbolicTangent::HyperbolicTangent()
{ }

HyperbolicTangent::~HyperbolicTangent()
{ }

double HyperbolicTangent::eval(double x) const
{
    if (m_disabled)
        return 0.0;

    return std::tanh(m_K * x);
}

double HyperbolicTangent::evalDerivative(double x) const
{
    if (m_disabled)
        return 0.0;

    double sechKx = 1.0/std::cosh(m_K * x);
    return m_K * sechKx * sechKx;
}

double HyperbolicTangent::evalDoubleDerivative(double x) const
{
    if (m_disabled)
        return 0.0;

    double sechKx = 1.0/std::cosh(m_K * x);
    return -2.0 * m_K * m_K * eval(x) * sechKx * sechKx;
}
