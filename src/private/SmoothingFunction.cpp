/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Utilities/SmoothingFunction.h>
#include <cmath>

using namespace DynamicalPlanner::Private;

SmoothingFunction::SmoothingFunction()
    : m_K(1.0)
      , m_disabled(false)
      , m_disabledValue(1.0)
{ }

SmoothingFunction::~SmoothingFunction()
{ }

void SmoothingFunction::setScaling(double k)
{
    m_K = k;
}

void SmoothingFunction::disable(double constantValue)
{
    m_disabled = true;
    m_disabledValue = constantValue;
}
