/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_ABSOLUTEVELOCITYEXPRESSION_H
#define DPLANNER_ABSOLUTEVELOCITYEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        levi::Expression AbsoluteLeftVelocityExpression(ExpressionsServer* expressionsServer,
                                                        const levi::Variable& baseTwist,
                                                        const std::string &targetFrame);

        levi::Expression AbsoluteLeftVelocityJointsDerivativeExpression(ExpressionsServer* expressionsServer,
                                                                        const levi::Variable& baseTwist,
                                                                        const std::string &targetFrame);

    }
}

#endif // DPLANNER_ABSOLUTEVELOCITYEXPRESSION_H
