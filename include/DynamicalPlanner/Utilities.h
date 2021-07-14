/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_UTILITIES_H
#define DPLANNER_UTILITIES_H

#include <DynamicalPlanner/Settings.h>
#include <DynamicalPlanner/RectangularFoot.h>
#include <DynamicalPlanner/State.h>

namespace DynamicalPlanner {
namespace Utilities {

bool FillDefaultInitialState(const DynamicalPlanner::Settings& inputSettings, const iDynTree::VectorDynSize& desiredJoints,
                             DynamicalPlanner::RectangularFoot &leftFoot, DynamicalPlanner::RectangularFoot &rightFoot,
                             DynamicalPlanner::State& initialState);

}
}

#endif // DPLANNER_UTILITIES_H
