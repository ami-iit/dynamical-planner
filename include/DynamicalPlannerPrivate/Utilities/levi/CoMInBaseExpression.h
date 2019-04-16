/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_COMPOSITIONEXPRESSION_H
#define DPLANNER_COMPOSITIONEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <DynamicalPlannerPrivate/Utilities/levi/TransformExpression.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        levi::Expression CoMInBaseExpression(ExpressionsServer* expressionsServer);

        levi::Expression CoMInBaseJointsDoubleDerivative(ExpressionsServer* expressionsServer, long column);

        levi::Expression CoMAdjointTransformWrench(const levi::Expression& worldToBaseRotation, const levi::Expression& comInBasePosition);
    }
}

#endif // DPLANNER_COMPOSITIONEXPRESSION_H
