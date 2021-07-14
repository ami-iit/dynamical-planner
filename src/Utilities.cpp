/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Utilities.h>
#include <iDynTree/KinDynComputations.h>
#include <cmath>

namespace DynamicalPlanner {
namespace Utilities {

bool FillDefaultInitialState(const Settings &inputSettings, const iDynTree::VectorDynSize &desiredJoints, RectangularFoot &leftFoot, RectangularFoot &rightFoot, State &initialState)
{
    if (!inputSettings.isValid())
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::FillDefaultInitialState] The input settings are not valid." <<std::endl;
        return false;
    }

    const DynamicalPlanner::SettingsStruct& settings = inputSettings.getSettings();

    iDynTree::KinDynComputations kinDyn;

    if (!kinDyn.loadRobotModel(settings.robotModel))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::FillDefaultInitialState] Failed to load the robot model in the KinDynComputations object." <<std::endl;
        return false;
    }

    initialState.resize(desiredJoints.size(), settings.leftPointsPosition.size());

    const iDynTree::Model& model = kinDyn.model();

    if (!kinDyn.setFloatingBase(model.getLinkName(model.getFrameLink(model.getFrameIndex(settings.leftFrameName)))))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::FillDefaultInitialState] Failed to set the floating base to the KinDynComputations object." <<std::endl;
        return false;
    }

    iDynTree::Vector3 gravity;
    gravity.zero();
    gravity(2) = -9.81;

    iDynTree::Twist baseVelocity = iDynTree::Twist::Zero();

    if (!kinDyn.setRobotState(model.getFrameTransform(model.getFrameIndex(settings.leftFrameName)).inverse(), desiredJoints,
                              baseVelocity, iDynTree::VectorDynSize(desiredJoints.size()), gravity))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::FillDefaultInitialState] Failed to set the robot state in the KinDynComputations object." <<std::endl;
        return false;
    }

    initialState.comPosition = kinDyn.getCenterOfMassPosition();

    initialState.jointsConfiguration = desiredJoints;

    kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
    initialState.momentumInCoM = kinDyn.getCentroidalTotalMomentum().asVector();

    initialState.time = 0.0;

    initialState.worldToBaseTransform = kinDyn.getWorldTransform(settings.floatingBaseName);

    iDynTree::Transform leftTransform = kinDyn.getWorldTransform(settings.leftFrameName);
    iDynTree::Transform rightTransform = kinDyn.getWorldTransform(settings.rightFrameName);

    double totalMass = 0.0;

    for(size_t l=0; l < model.getNrOfLinks(); l++)
    {
        totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
    }

    double normalForce = totalMass * 9.81;

    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointPosition = leftTransform * settings.leftPointsPosition[i];
        initialState.rightContactPointsState[i].pointPosition = rightTransform * settings.rightPointsPosition[i];
    }

    iDynTree::Wrench leftWrench, rightWrench;
    leftWrench.zero();
    rightWrench.zero();

    leftWrench(2) = normalForce/(1 + std::fabs(initialState.comPosition(1) - leftTransform.getPosition()(1))/std::fabs(initialState.comPosition(1) - rightTransform.getPosition()(1)));

    rightWrench(2) = normalForce - leftWrench(2);

    leftWrench(4) = -leftWrench(2) * (initialState.comPosition(0) - leftTransform.getPosition()(0));
    rightWrench(4) = -rightWrench(2) * (initialState.comPosition(0) - rightTransform.getPosition()(0));

    std::vector<iDynTree::Force> leftPointForces, rightPointForces;

    if (!leftFoot.getForces(leftWrench, leftPointForces))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::FillDefaultInitialState] Failed to compute the forces for the left foot." <<std::endl;
        return false;
    }

    if (!rightFoot.getForces(rightWrench, rightPointForces))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::FillDefaultInitialState] Failed to compute the forces for the right foot." <<std::endl;
        return false;
    }

    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointForce = leftPointForces[i];
        initialState.rightContactPointsState[i].pointForce =rightPointForces[i];
    }

    return true;
}

}
}
