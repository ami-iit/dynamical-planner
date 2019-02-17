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
#include <DynamicalPlannerPrivate/Utilities/levi/TransformExpression.h>
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

    levi::Variable* baseQuaternionVelocity();

    levi::Expression* baseTwist();

    levi::Variable* jointsPosition();

    levi::Variable* jointsVelocity();

    TransformExpression* worldToBase();

    levi::Expression* comInBase();

    levi::Expression* adjointTransform(const std::string& baseFrame,
                                       const std::string &targetFrame);

    levi::Expression* adjointTransformWrench(const std::string& baseFrame,
                                             const std::string &targetFrame);

    levi::Expression* relativePosition(const std::string& baseFrame,
                                       const std::string &targetFrame);

    TransformExpression* relativeTransform(const std::string& baseFrame,
                                           const std::string &targetFrame);

    levi::Expression* relativeVelocity(const std::string& baseFrame,
                                       const std::string &targetFrame);

};

#endif // DPLANNER_EXPRESSIONSSERVER_H
