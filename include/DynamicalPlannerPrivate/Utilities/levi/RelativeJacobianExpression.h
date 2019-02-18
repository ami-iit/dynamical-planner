/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_RELATIVEJACOBIANEXPRESSION_H
#define DPLANNER_RELATIVEJACOBIANEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        levi::Expression RelativeLeftJacobianExpression(ExpressionsServer* expressionsServer,
                                                        const std::string& baseFrame,
                                                        const std::string &targetFrame);

    }
}

#endif // DPLANNER_RELATIVEJACOBIANEXPRESSION_H
