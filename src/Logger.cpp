/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Logger.h>
#include <matlogger2/matlogger2.h>
#include <iDynTree/Core/EigenHelpers.h>

using namespace DynamicalPlanner;

template <typename Vector>
void addToLogger(XBot::MatLogger2::Ptr logger, const std::string& variableName, const Vector& value) {
    logger->add(variableName, iDynTree::toEigen(value));
}

void Logger::saveSolutionVectorsToFile(const std::string &matFileName, const std::vector<State> &states, const std::vector<Control> &controls)
{
    auto logger = XBot::MatLogger2::MakeLogger(matFileName);

    int statesSize = static_cast<int>(states.size());
    int controlsSize = static_cast<int>(controls.size());

    if (statesSize) {
        for (size_t point = 0; point < states.front().leftContactPointsState.size(); ++point) {
            logger->create("leftPoint" + std::to_string(point) + "Force", 3, 1, statesSize);
            logger->create("leftPoint" + std::to_string(point) + "Position", 3, 1, statesSize);
        }

        for (size_t point = 0; point < states.front().rightContactPointsState.size(); ++point) {
            logger->create("rightPoint" + std::to_string(point) + "Force", 3, 1, statesSize);
            logger->create("rightPoint" + std::to_string(point) + "Position", 3, 1, statesSize);
        }

        logger->create("momentumInCoM", 6, 1, statesSize);
        logger->create("comPosition", 3, 1, statesSize);
        logger->create("basePosition", 3, 1, statesSize);
        logger->create("baseQuaternion", 4, 1, statesSize);
        logger->create("jointsConfiguraion", static_cast<int>(states.front().jointsConfiguration.size()), 1, statesSize);
        logger->create("stateTime", 1, 1, statesSize);
    }

    if (controlsSize) {
        for (size_t point = 0; point < controls.front().leftContactPointsControl.size(); ++point) {
            logger->create("leftPoint" + std::to_string(point) + "ForceControl", 3, 1, controlsSize);
            logger->create("leftPoint" + std::to_string(point) + "VelocityControl", 3, 1, controlsSize);
        }

        for (size_t point = 0; point < controls.front().rightContactPointsControl.size(); ++point) {
            logger->create("rightPoint" + std::to_string(point) + "ForceControl", 3, 1, controlsSize);
            logger->create("rightPoint" + std::to_string(point) + "VelocityControl", 3, 1, controlsSize);
        }

        logger->create("baseLinearVelocity", 3, 1, controlsSize);
        logger->create("baseQuaternionDerivative", 4, 1, controlsSize);
        logger->create("jointsVelocity", static_cast<int>(controls.front().jointsVelocity.size()), 1, controlsSize);
        logger->create("controlTime", 1, 1, controlsSize);
    }

    for (auto state = states.begin(); state != states.end(); ++state) {
        for (size_t point = 0; point < states.front().leftContactPointsState.size(); ++point) {
            addToLogger(logger, "leftPoint" + std::to_string(point) + "Force", state->leftContactPointsState[point].pointForce);
            addToLogger(logger,"leftPoint" + std::to_string(point) + "Position", state->leftContactPointsState[point].pointPosition);
        }

        for (size_t point = 0; point < states.front().rightContactPointsState.size(); ++point) {
            addToLogger(logger, "rightPoint" + std::to_string(point) + "Force", state->rightContactPointsState[point].pointForce);
            addToLogger(logger,"rightPoint" + std::to_string(point) + "Position", state->rightContactPointsState[point].pointPosition);
        }

        addToLogger(logger, "momentumInCoM", state->momentumInCoM);
        addToLogger(logger, "comPosition", state->comPosition);
        addToLogger(logger, "basePosition", state->worldToBaseTransform.getPosition());
        addToLogger(logger, "baseQuaternion", state->worldToBaseTransform.getRotation().asQuaternion());
        addToLogger(logger, "jointsConfiguraion", state->jointsConfiguration);
        logger->add("stateTime", state->time);

        logger->flush_available_data();
    }

    for (auto control = controls.begin(); control != controls.end(); ++control) {
        for (size_t point = 0; point < controls.front().leftContactPointsControl.size(); ++point) {
            addToLogger(logger, "leftPoint" + std::to_string(point) + "ForceControl", control->leftContactPointsControl[point].pointForceControl);
            addToLogger(logger, "leftPoint" + std::to_string(point) + "VelocityControl",
                        control->leftContactPointsControl[point].pointVelocityControl);
        }

        for (size_t point = 0; point < controls.front().rightContactPointsControl.size(); ++point) {
            addToLogger(logger, "rightPoint" + std::to_string(point) + "ForceControl", control->rightContactPointsControl[point].pointForceControl);
            addToLogger(logger, "rightPoint" + std::to_string(point) + "VelocityControl",
                        control->rightContactPointsControl[point].pointVelocityControl);
        }

        addToLogger(logger, "baseLinearVelocity", control->baseLinearVelocity);
        addToLogger(logger, "baseQuaternionDerivative", control->baseQuaternionDerivative);
        addToLogger(logger, "jointsVelocity", control->jointsVelocity);
        logger->add("controlTime", control->time);

        logger->flush_available_data();
    }

    logger.reset(); // manually destroy the logger in order to flush to disk
}

