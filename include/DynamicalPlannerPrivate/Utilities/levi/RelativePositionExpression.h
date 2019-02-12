/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_RELATIVEPOSITIONEXPRESSION_H
#define DPLANNER_RELATIVEPOSITIONEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        levi::Expression RelativePositionExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                    RobotState *robotState,
                                                    const std::string& baseFrame,
                                                    const std::string &targetFrame,
                                                    levi::Variable jointsVariable,
                                                    levi::ScalarVariable timeVariable);

        levi::Expression RelativePositionExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                    RobotState *robotState,
                                                    const std::string& baseFrame,
                                                    const std::string &targetFrame,
                                                    levi::Variable jointsVariable,
                                                    levi::ScalarVariable timeVariable,
                                                    const levi::Expression& relativeJacobian,
                                                    const levi::Expression& relativeRotation);

    }
}

#endif // DPLANNER_RELATIVEPOSITIONEXPRESSION_H
