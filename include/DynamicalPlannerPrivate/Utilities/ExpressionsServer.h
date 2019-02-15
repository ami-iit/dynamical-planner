/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_EXPRESSIONSSERVER_H
#define DPLANNER_EXPRESSIONSSERVER_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class ExpressionsServer;
    }
}

class DynamicalPlanner::Private::ExpressionsServer {
    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    ExpressionsServer(std::shared_ptr<DynamicalPlanner::Private::TimelySharedKinDynComputations> timelySharedKinDyn);

    ~ExpressionsServer();

    bool updateRobotState(double time, const RobotState &currentState);

    levi::Expression* baseRotation();

    levi::Variable* baseQuaternion(); //not normalized

    levi::Variable* basePosition();

    levi::Variable* baseLinearVelocity();

    levi::Variable* baseAngularVelocity();

    levi::Variable* jointsPosition();

    levi::Variable* jointsVelocity();

    levi::Expression* adjointTransform(const std::string& baseFrame,
                                       const std::string &targetFrame);

    levi::Expression* adjointTransformWrench(const std::string& baseFrame,
                                             const std::string &targetFrame);

};

#endif // DPLANNER_EXPRESSIONSSERVER_H
