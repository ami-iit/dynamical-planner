/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_COMPOSITIONEXPRESSION_H
#define DPLANNER_COMPOSITIONEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Utilities/levi/TransformExpression.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

//        levi::Expression CoMPositionExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
//                                               RobotState *robotState,
//                                               const TransformExpression& baseTransform,
//                                               levi::Variable jointsVariable,
//                                               levi::ScalarVariable timeVariable);

        levi::Expression CoMMixedJacobianExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                    RobotState *robotState,
                                                    const TransformExpression& baseTransform,
                                                    levi::Variable jointsVariable,
                                                    levi::ScalarVariable timeVariable);

    }
}

#endif // DPLANNER_COMPOSITIONEXPRESSION_H
