/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_VISUALIZER_H
#define DPLANNER_VISUALIZER_H

#include <DynamicalPlanner/State.h>
#include <iDynTree/Model/Model.h>
#include <iDynTree/Core/Direction.h>
#include <iDynTree/Core/Position.h>
#include <memory>
#include <functional>

namespace DynamicalPlanner {
    class Visualizer;
}

class DynamicalPlanner::Visualizer {
public:

    Visualizer();

    ~Visualizer();

    bool setModel(const iDynTree::Model& model);

    bool visualizeState(const DynamicalPlanner::State& stateToVisualize);

    bool visualizeStates(const std::vector<State>& states, double endTime = -1.0);

    bool visualizeStates(const std::vector<State>& states,
                         const std::vector<iDynTree::Position>& cameraPosition,
                         const std::vector<iDynTree::Position>& cameraTarget,
                         double endTime = -1.0);

    bool visualizeStatesAndSaveAnimation(const std::vector<State>& states, const std::string& workingFolder,
                                const std::string& fileName, const std::string& fileExtension = "gif", double endTime = -1.0);

    bool setCameraPosition(const iDynTree::Position& cameraPosition);

    bool setCameraTarget(const iDynTree::Position& cameraTarget);

    bool setLightDirection(const iDynTree::Direction& lightDirection);

    void visualizeWorldFrame(bool visualizeWorldFrame = true);

    bool visualizeMPCStatesAndSaveAnimation(const std::vector<State>& states, std::function<iDynTree::Position(const State&)> stateCameraControl,
                                            const std::vector<std::vector<State>>& fullStates, std::function<iDynTree::Position(const State&)> fullStateCameraControl,
                                            const std::string& workingFolder, const std::string& fileName, const std::string& fileExtension = "gif",
                                            double speedUp = 1.0, double endTime = -1.0);

private:

    class VisualizerImplementation;
    std::unique_ptr<VisualizerImplementation> m_pimpl;
};

#endif // DPLANNER_VISUALIZER_H
