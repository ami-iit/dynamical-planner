/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_COSTS_H
#define DPLANNER_COSTS_H

#include <DynamicalPlannerPrivate/Costs/ForceMeanCost.h>
#include <DynamicalPlannerPrivate/Costs/FrameOrientationCost.h>
#include <DynamicalPlannerPrivate/Costs/StaticTorquesCost.h>
#include <DynamicalPlannerPrivate/Costs/SwingCost.h>
#include <DynamicalPlannerPrivate/Costs/PhantomForcesCost.h>
#include <DynamicalPlannerPrivate/Costs/MeanPointPositionCost.h>
#include <DynamicalPlannerPrivate/Costs/CoMVelocityCost.h>
#include <DynamicalPlannerPrivate/Costs/FootYawCost.h>
#include <DynamicalPlannerPrivate/Costs/FeetDistanceCost.h>
#include <DynamicalPlannerPrivate/Costs/JointsVelocityForPosturalCost.h>
#include <DynamicalPlannerPrivate/Costs/ComplementarityCost.h>
#include <iDynTree/L2NormCost.h>

#endif // DPLANNER_COSTS_H
