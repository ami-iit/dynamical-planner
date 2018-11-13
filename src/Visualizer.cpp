/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Visualizer.h>
#include <iDynTree/Visualizer.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>


using namespace DynamicalPlanner;

class Visualizer::VisualizerImplementation {
public:

    iDynTree::Visualizer viz;
    iDynTree::Position defaultCameraPosition, defaultCameraTarget;

    VisualizerImplementation() {}
    ~VisualizerImplementation(){}
};

Visualizer::Visualizer()
    : m_pimpl(std::make_unique<VisualizerImplementation>())
{
    m_pimpl->viz.init();
    setCameraPosition(iDynTree::Position(1.0, 0.0, 0.5));
    setCameraTarget(iDynTree::Position(0.4, 0.0, 0.5));
    double sqrt2 = std::sqrt(2.0);
    setLightDirection(iDynTree::Direction(-0.5/sqrt2, 0, -0.5/sqrt2));
    m_pimpl->viz.vectors().setVectorsAspect(0.01, 0.0, 0.01);
}

Visualizer::~Visualizer()
{
    m_pimpl->viz.close();
}

bool Visualizer::setModel(const iDynTree::Model &model)
{
    if (m_pimpl->viz.getNrOfVisualizedModels()) {
        std::cerr << "[ERROR][Visualizer::setModel] Model already set." << std::endl;
        return false;
    }
    bool modelLoaded = m_pimpl->viz.addModel(model, "DynamicalPlannerVisualizer");

    if (!modelLoaded) {
        std::cerr << "[ERROR][Visualizer::setModel] Failed to set model for visualization." << std::endl;
        return false;
    }

    return modelLoaded;
}

bool Visualizer::visualizeState(const State &stateToVisualize)
{
    if (!(m_pimpl->viz.getNrOfVisualizedModels())) {
        std::cerr << "[ERROR][Visualizer::visualizeState] First you have to load a model." << std::endl;
        return false;
    }

    m_pimpl->viz.modelViz(0).setPositions(stateToVisualize.worldToBaseTransform, stateToVisualize.jointsConfiguration);

    iDynTree::IVectorsVisualization& forcesViz = m_pimpl->viz.vectors();

    size_t vectorIndex = 0;
    iDynTree::Position posBuf;

    for (const ContactPointState& point : stateToVisualize.leftContactPointsState) {
        iDynTree::toEigen(posBuf) = iDynTree::toEigen(point.pointPosition);
        if (vectorIndex < forcesViz.getNrOfVectors()) {
            forcesViz.updateVector(vectorIndex, posBuf, point.pointForce);
        } else {
            forcesViz.addVector(posBuf, point.pointForce);
        }
        vectorIndex++;
    }

    for (const ContactPointState& point : stateToVisualize.rightContactPointsState) {
        iDynTree::toEigen(posBuf) = iDynTree::toEigen(point.pointPosition);
        if (vectorIndex < forcesViz.getNrOfVectors()) {
            forcesViz.updateVector(vectorIndex, posBuf, point.pointForce);
        } else {
            forcesViz.addVector(posBuf, point.pointForce);
        }
        vectorIndex++;
    }

    m_pimpl->viz.draw();

    return true;
}

bool Visualizer::visualizeStates(const std::vector<State> &states, double endTime)
{
    if (!(m_pimpl->viz.getNrOfVisualizedModels())) {
        std::cerr << "[ERROR][Visualizer::visualizeState] First you have to load a model." << std::endl;
        return false;
    }

    for (size_t i = 0; i < states.size(); ++i) {
        if ((endTime > 0) && (states[i].time > endTime)) {
            break;
        }
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        visualizeState(states[i]);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if ((i + 1) < states.size()) {
            std::chrono::milliseconds durationMs(static_cast<int>(std::round((states[i + 1].time - states[i].time)*1000)));
            std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

            if (elapsed < durationMs) {
                std::this_thread::sleep_for (durationMs - elapsed);
            }
        }
    }

    return true;
}

bool Visualizer::visualizeStates(const std::vector<State> &states, const std::vector<iDynTree::Position> &cameraPosition, const std::vector<iDynTree::Position> &cameraTarget, double endTime)
{
    if (!(m_pimpl->viz.getNrOfVisualizedModels())) {
        std::cerr << "[ERROR][Visualizer::visualizeState] First you have to load a model." << std::endl;
        return false;
    }

    if (cameraPosition.size() != states.size()) {
        std::cerr << "[ERROR][Visualizer::visualizeState] The cameraPosition vector should have the same dimension of the states vector." << std::endl;
        return false;
    }

    if (cameraTarget.size() != states.size()) {
        std::cerr << "[ERROR][Visualizer::visualizeState] The cameraTarget vector should have the same dimension of the states vector." << std::endl;
        return false;
    }

    for (size_t i = 0; i < states.size(); ++i) {
        if ((endTime > 0) && (states[i].time > endTime)) {
            break;
        }
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        m_pimpl->viz.camera().setPosition(cameraPosition[i]);
        m_pimpl->viz.camera().setTarget(cameraTarget[i]);
        visualizeState(states[i]);
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if ((i + 1) < states.size()) {
            std::chrono::milliseconds durationMs(static_cast<int>(std::round((states[i + 1].time - states[i].time)*1000)));
            std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

            if (elapsed < durationMs) {
                std::this_thread::sleep_for (durationMs - elapsed);
            }
        }
    }

    m_pimpl->viz.camera().setPosition(m_pimpl->defaultCameraPosition);

    m_pimpl->viz.camera().setTarget(m_pimpl->defaultCameraTarget);

    return true;
}

bool Visualizer::visualizeStatesAndSaveAnimation(const std::vector<State> &states, const std::string &workingFolder, const std::string &fileName, const std::string &fileExtension, double endTime)
{
    if (!(m_pimpl->viz.getNrOfVisualizedModels())) {
        std::cerr << "[ERROR][Visualizer::visualizeState] First you have to load a model." << std::endl;
        return false;
    }

    unsigned int digits = static_cast<unsigned int>(std::floor(std::log(states.size()) + 1));

    size_t i = 0;

    while (i < states.size() && (!(endTime < 0) || (states[i].time <= endTime))) {

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        visualizeState(states[i]);
        m_pimpl->viz.drawToFile(workingFolder + "/" + fileName + "_img_" + std::string(digits - std::to_string(i).size(), '0') + std::to_string(i) + ".png");
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if ((i + 1) < states.size()) {
            std::chrono::milliseconds durationMs(static_cast<int>(std::round((states[i + 1].time - states[i].time)*1000)));
            std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

            if (elapsed < durationMs) {
                std::this_thread::sleep_for (durationMs - elapsed);
            }
        }

        ++i;
    }

    if (i == 0) {
        return true;
    }

    std::string fps = std::to_string(static_cast<int>(std::round(i/states[i-1].time)));

    auto frameArgs          = " -framerate " + fps;
    auto inputArgs          = " -i " + workingFolder + "/" + fileName + "_img_%0" + std::to_string(digits) + "d" + ".png";
    auto overwriteArgs      = " -y";
    auto outputArgs         = " " + workingFolder + "/" + fileName + "." + fileExtension;
    auto pixFormatArgs      = "";

    // http://superuser.com/questions/533695/how-can-i-convert-a-series-of-png-images-to-a-video-for-youtube#answers-header
    if (fileExtension == "mp4") {
        pixFormatArgs = " -pix_fmt yuv420p";
    }

    auto args = "ffmpeg" + frameArgs + inputArgs + overwriteArgs + pixFormatArgs + outputArgs;

    std::cout << "[INFO][Visualizer::visualizeState] Generating video with the following arguments:\n" << args << std::endl;

    int useless_int = system(args.c_str());

    return true;
}

bool Visualizer::setCameraPosition(const iDynTree::Position &cameraPosition)
{
    m_pimpl->viz.camera().setPosition(cameraPosition);
    m_pimpl->defaultCameraPosition = cameraPosition;
    return true;
}

bool Visualizer::setCameraTarget(const iDynTree::Position &cameraTarget)
{
    m_pimpl->viz.camera().setTarget(cameraTarget);
    m_pimpl->defaultCameraTarget = cameraTarget;
    return true;
}

bool Visualizer::setLightDirection(const iDynTree::Direction &lightDirection)
{
    m_pimpl->viz.enviroment().lightViz("sun").setDirection(lightDirection);

    return true;

}
