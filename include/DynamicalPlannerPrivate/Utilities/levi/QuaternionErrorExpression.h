/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_QUATERNIONERROREXPRESSION_H
#define DPLANNER_QUATERNIONERROREXPRESSION_H

#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <levi/ForwardDeclarations.h>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        levi::Expression QuaternionError(const std::string &desiredFrame, ExpressionsServer* expressionsServer, const levi::Variable& desiredQuaternion);
    }
}

#endif // DPLANNER_QUATERNIONERROREXPRESSION_H
