/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_ADJOINTTRANSFORMEXPRESSION_H
#define DPLANNER_ADJOINTTRANSFORMEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        levi::Expression AdjointTransformExpression(ExpressionsServer* expressionServer,
                                                    const std::string& baseFrame,
                                                    const std::string &targetFrame);

        levi::Expression AdjointTransformExpressionJointsDerivative(ExpressionsServer* expressionsServer,
                                                                    const std::string& baseFrame,
                                                                    const std::string &targetFrame,
                                                                    long column);

        levi::Expression AdjointTransformWrenchExpression(ExpressionsServer* expressionServer,
                                                          const std::string& baseFrame,
                                                          const std::string &targetFrame);

        levi::Expression AdjointTransformWrenchExpressionJointsDerivative(ExpressionsServer* expressionServer,
                                                                          const std::string& baseFrame,
                                                                          const std::string &targetFrame,
                                                                          long column);

    }
}

#endif // DPLANNER_ADJOINTTRANSFORMEXPRESSION_H
