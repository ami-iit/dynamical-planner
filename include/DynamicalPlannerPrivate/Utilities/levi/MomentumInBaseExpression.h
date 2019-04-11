/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_MOMENTUMINBASEEXPRESSION_H
#define DPLANNER_MOMENTUMINBASEEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        levi::Expression MomentumInBaseExpression(ExpressionsServer* expressionsServer,
                                                  const levi::Variable& baseTwist);

        levi::Expression MomentumInBaseExpressionJointsDerivativeExpression(ExpressionsServer* expressionsServer,
                                                                            const levi::Variable& baseTwist);

        levi::Expression MomentumInBaseExpressionBaseTwistDerivativeExpression(ExpressionsServer* expressionsServer,
                                                                               const levi::Expression& baseTwistJacobian);

        levi::Expression MomentumInBaseJointsVelocityDerivativeExpression(ExpressionsServer* expressionsServer);

        levi::Expression MomentumInBaseBaseTwistJointsDerivativeExpression(ExpressionsServer* expressionsServer,
                                                                           const levi::Expression& baseTwistJacobian,
                                                                           long column);

        levi::Expression MomentumInBaseExpressionJointsDoubleDerivativeExpression(ExpressionsServer* expressionsServer,
                                                                                  const levi::Variable& baseTwist,
                                                                                  long column);

    }
}

#endif // DPLANNER_MOMENTUMINBASEEXPRESSION_H
