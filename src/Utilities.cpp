/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Utilities.h>
#include <iDynTree/KinDynComputations.h>
#include <iDynTree/Core/EigenHelpers.h>
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

bool SetMinContactPointToZero_impl(iDynTree::KinDynComputations& kinDyn, const DynamicalPlanner::SettingsStruct& settings, DynamicalPlanner::State &initialState) {

    iDynTree::Vector3 gravity;
    gravity.zero();
    gravity(2) = -9.81;

    if (!kinDyn.setRobotState(initialState.worldToBaseTransform, initialState.jointsConfiguration,
                              iDynTree::Twist::Zero(), iDynTree::VectorDynSize(initialState.jointsConfiguration.size()), gravity))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::SetMinContactPointToZero] Failed to set the robot state in the KinDynComputations object." <<std::endl;
        return false;
    }

    initialState.comPosition = kinDyn.getCenterOfMassPosition();

    iDynTree::Transform leftTransform = kinDyn.getWorldTransform(settings.leftFrameName);
    iDynTree::Transform rightTransform = kinDyn.getWorldTransform(settings.rightFrameName);

    double minZ = 1.0;

    for (size_t i = 0; i < settings.leftPointsPosition.size(); ++i) {
        initialState.leftContactPointsState[i].pointPosition = leftTransform * settings.leftPointsPosition[i];

        if (initialState.leftContactPointsState[i].pointPosition(2) < minZ) {
            minZ = initialState.leftContactPointsState[i].pointPosition(2);
        }

        initialState.rightContactPointsState[i].pointPosition = rightTransform * settings.rightPointsPosition[i];

        if (initialState.rightContactPointsState[i].pointPosition(2) < minZ) {
            minZ = initialState.rightContactPointsState[i].pointPosition(2);
        }
    }

    if (!iDynTree::checkDoublesAreEqual(minZ, 0.0, 1e-6)) {
        iDynTree::Position initialPosition = initialState.worldToBaseTransform.getPosition();
        initialPosition(2) -= minZ;
        initialState.worldToBaseTransform.setPosition(initialPosition);
        return SetMinContactPointToZero_impl(kinDyn, settings, initialState);
    }

    return true;
}

bool SetMinContactPointToZero(const Settings &inputSettings, State &initialState) {

    if (!inputSettings.isValid())
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::SetMinContactPointToZero] The input settings are not valid." <<std::endl;
        return false;
    }

    const DynamicalPlanner::SettingsStruct& settings = inputSettings.getSettings();

    iDynTree::KinDynComputations kinDyn;

    if (!kinDyn.loadRobotModel(settings.robotModel))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::SetMinContactPointToZero] Failed to load the robot model in the KinDynComputations object." <<std::endl;
        return false;
    }

    if (!kinDyn.setFloatingBase(settings.floatingBaseName))
    {
        std::cerr << "[ERROR][DynamicalPlanner::Utilities::SetMinContactPointToZero] Failed to set the floating base to the KinDynComputations object." <<std::endl;
        return false;
    }

    return SetMinContactPointToZero_impl(kinDyn, settings, initialState);
}

TranslatingCoMStateGuess::TranslatingCoMStateGuess(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> comReference, const State &initialState)
    : m_state(initialState)
    , m_initialState(initialState)
    , m_comReference(comReference)
{ }

TranslatingCoMStateGuess::~TranslatingCoMStateGuess()
{ }

State &TranslatingCoMStateGuess::get(double time, bool &isValid) {
    iDynTree::toEigen(m_state.comPosition) = iDynTree::toEigen(m_comReference->get(time, isValid));
    m_state.jointsConfiguration = m_initialState.jointsConfiguration;
    m_state.momentumInCoM.zero();
    m_state.worldToBaseTransform.setRotation(m_initialState.worldToBaseTransform.getRotation());
    iDynTree::Position basePosition, comDifference;
    iDynTree::toEigen(comDifference) = iDynTree::toEigen(m_state.comPosition) - iDynTree::toEigen(m_initialState.comPosition);
    iDynTree::toEigen(basePosition) = iDynTree::toEigen(m_initialState.worldToBaseTransform.getPosition()) + iDynTree::toEigen(comDifference);
    m_state.worldToBaseTransform.setPosition(basePosition);

    for (size_t i = 0; i < m_state.leftContactPointsState.size(); ++i) {
        iDynTree::toEigen(m_state.leftContactPointsState[i].pointPosition) = iDynTree::toEigen(m_initialState.leftContactPointsState[i].pointPosition) + iDynTree::toEigen(comDifference);
    }

    for (size_t i = 0; i < m_state.rightContactPointsState.size(); ++i) {
        iDynTree::toEigen(m_state.rightContactPointsState[i].pointPosition) = iDynTree::toEigen(m_initialState.rightContactPointsState[i].pointPosition) + iDynTree::toEigen(comDifference);
    }

    isValid = true;

    m_state.time = time;
    return m_state;
}

bool SimpleWalkingStateMachine::initialize(const iDynTree::Vector3 startingReference, const iDynTree::Vector3 &stepIncrement,
                                           double stepDuration, double horizon, double minimumDt,
                                           double weightIncreaseX, double weightIncreaseY, double weightIncreaseZ,
                                           double forceThreshold)
{
    if (forceThreshold < 0)
    {
        std::cerr << "[ERROR][SimpleWalkingStateMachine::initialize] The force threshold is supposed to be non-negative." << std::endl;
        return false;
    }

    if (minimumDt > horizon || minimumDt > stepDuration)
    {
        std::cerr << "[ERROR][SimpleWalkingStateMachine::initialize] The minimumDt is supposed to be greater than the horizon and the step duration." << std::endl;
        return false;
    }

    m_references = std::make_shared<PositionReferenceGenerator>(2, weightIncreaseX, weightIncreaseY, weightIncreaseZ);
    iDynTree::toEigen(m_stepIncrement) = iDynTree::toEigen(stepIncrement);
    m_stepDuration = stepDuration;
    m_horizon = horizon;
    m_minimumDt = minimumDt;
    m_forceThreshold = forceThreshold;

    iDynTree::toEigen(m_references->at(0).desiredPosition) = iDynTree::toEigen(startingReference);
    m_references->at(0).activeRange.setTimeInterval(0.0, stepDuration);
    m_references->at(1).desiredPosition = m_references->at(0).desiredPosition;
    m_references->at(1).activeRange.setTimeInterval(m_horizon + 1.0, m_horizon + 1.0);

    return true;
}

void SimpleWalkingStateMachine::setVerbose(bool verbose)
{
    m_verbose = verbose;
}

bool SimpleWalkingStateMachine::advance(const State &currentState, const State &futureState)
{
    if (!m_references)
    {
        std::cerr << "[ERROR][SimpleWalkingStateMachine::advance] The initialize method was either not called or did not return correctly." << std::endl;
        return false;
    }

    double futureMeanPositionError = (iDynTree::toEigen(futureState.computeFeetCentroid()) - iDynTree::toEigen(m_references->at(0).desiredPosition)).norm();

    if (m_verbose)
    {
        std::cout << "[SimpleWalkingStateMachine::advance] Future mean point error: " << futureMeanPositionError << std::endl;
    }

    double meanPositionError = (iDynTree::toEigen(currentState.computeFeetCentroid()) - iDynTree::toEigen(m_references->at(0).desiredPosition)).norm();
    if (m_verbose)
    {
        std::cout << "[SimpleWalkingStateMachine::advance] Mean point error: " << meanPositionError << std::endl;
    }

    double minimumForce = currentState.minimumPointForceOnForwardFoot(m_stepIncrement);
    if (m_verbose)
    {
        std::cout << "[SimpleWalkingStateMachine::advance] Minimum force: " << minimumForce << std::endl;
    }

    if ((futureMeanPositionError < 5e-3) && (m_references->at(1).activeRange.initTime() > m_horizon) && (m_references->at(0).activeRange.endTime() < (0.6 * m_horizon))) {
        //You are here if the step is already completed in the future. The end time of the current step has to be early enough to insert the new step at the end of the horizon.
        m_references->at(1).activeRange.setTimeInterval(m_horizon, m_horizon + m_stepDuration);
        m_references->at(1).desiredPosition = m_references->at(0).desiredPosition + m_stepIncrement;
        if (m_verbose)
        {
            std::cout << "[SimpleWalkingStateMachine::advance] Setting new position (" << m_references->at(1).desiredPosition.toString() << ") at the end of the horizon." << std::endl;
        }
    }

    if (m_references->at(1).activeRange.initTime() <= m_horizon) {

        double stepStart = m_references->at(0).activeRange.initTime() - currentState.time;

        if ((stepStart < 0.3*m_horizon) && (minimumForce < m_forceThreshold)) {
        // You are here if the second step is about to substitute the first, but the force in the forward foot is still too low
            m_references->at(0).activeRange.setTimeInterval(stepStart, m_references->at(0).activeRange.endTime());
            if (m_verbose)
            {
                std::cout << "[SimpleWalkingStateMachine::advance] New first step interval: [" << m_references->at(0).activeRange.initTime() << ", " << m_references->at(0).activeRange.endTime() << "]." << std::endl;
                std::cout << "[SimpleWalkingStateMachine::advance] Second step is kept constant: [" << m_references->at(1).activeRange.initTime() << ", " << m_references->at(1).activeRange.endTime() << "]." << std::endl;
            }

        } else if ((stepStart +  m_stepDuration) < m_minimumDt) {
            // You are here when the robot is ready to perform a new step.
            m_references->at(0).desiredPosition = m_references->at(1).desiredPosition;
            m_references->at(0).activeRange.setTimeInterval(0.0, m_stepDuration);
            m_references->at(1).activeRange.setTimeInterval(m_horizon + 1.0, m_horizon + 1.0);
            if (m_verbose)
            {
                std::cout << "[SimpleWalkingStateMachine::advance] Second step is now first step." << std::endl;
            }
        } else {
            // You are here if the force on the forward foot is high enough, but it is too early to perform a new step
            m_references->at(0).activeRange.setTimeInterval(stepStart, stepStart + m_stepDuration);
            std::cout << "New first step interval: [" << m_references->at(0).activeRange.initTime() << ", " << m_references->at(0).activeRange.endTime() << "]." << std::endl;
            stepStart = m_references->at(1).activeRange.initTime() - currentState.time;
            m_references->at(1).activeRange.setTimeInterval(stepStart, stepStart + m_stepDuration);
            if (m_verbose)
            {
                std::cout << "[SimpleWalkingStateMachine::advance] New second step interval: [" << m_references->at(1).activeRange.initTime() << ", " << m_references->at(1).activeRange.endTime() << "]." << std::endl;
            }
        }
    } else {
        // You are here if you are performing a step and it is not done yet. Hence the second step is not planned yet

        double stepStart = m_references->at(0).activeRange.initTime() - currentState.time;
        m_references->at(0).activeRange.setTimeInterval(stepStart, std::max(stepStart + m_stepDuration, 0.3 * m_horizon));
        if (m_verbose)
        {
            std::cout << "[SimpleWalkingStateMachine::advance] New first step interval: [" << m_references->at(0).activeRange.initTime() << ", " << m_references->at(0).activeRange.endTime() << "]." << std::endl;
        }
    }

    return true;
}

std::shared_ptr<PositionReferenceGenerator> SimpleWalkingStateMachine::references()
{
    return m_references;
}

}
}
