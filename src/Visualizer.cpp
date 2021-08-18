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
#include <iomanip>
#include <sstream>


using namespace DynamicalPlanner;

class Visualizer::VisualizerImplementation {
public:

    iDynTree::Visualizer viz;
    iDynTree::ITexture* textureForScreenshots;
    iDynTree::Position defaultCameraPosition, defaultCameraTarget;

    VisualizerImplementation() {}
    ~VisualizerImplementation(){}

    void updateViz(const State &stateToVisualize)
    {
        viz.modelViz(0).setPositions(stateToVisualize.worldToBaseTransform, stateToVisualize.jointsConfiguration);

        iDynTree::IVectorsVisualization& forcesViz = viz.vectors();

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
    }
};

Visualizer::Visualizer()
    : m_pimpl(std::make_unique<VisualizerImplementation>())
{
    iDynTree::VisualizerOptions options;
    options.winWidth = 1920;
    options.winHeight = 1080;
    m_pimpl->viz.init(options);
    setCameraPosition(iDynTree::Position(1.0, 0.0, 0.5));
    setCameraTarget(iDynTree::Position(0.4, 0.0, 0.5));
    double sqrt2 = std::sqrt(2.0);
    m_pimpl->viz.vectors().setVectorsAspect(0.01, 0.0, 0.01);
    m_pimpl->textureForScreenshots = m_pimpl->viz.textures().add("screenshotSaver", options);
    setLightDirection(iDynTree::Direction(-0.5/sqrt2, 0, -0.5/sqrt2));
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

    m_pimpl->updateViz(stateToVisualize);

    m_pimpl->viz.draw();
    m_pimpl->viz.run(); //This is to make sure that the window sizes are updated. This is done after the draw in case there are multiple visualizers. In fact, the draw method selects the correct window

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

    unsigned int digits = static_cast<unsigned int>(std::floor(std::log10(states.size()) + 1));

    size_t i = 0;

    while (i < states.size() && (!(endTime < 0) || (states[i].time <= endTime))) {

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        visualizeState(states[i]);
        //Using a texture as the dimension should be always the same and not depending on the window size
        m_pimpl->textureForScreenshots->drawToFile(workingFolder + "/" + fileName + "_img_" + std::string(digits - std::to_string(i).size(), '0') + std::to_string(i) + ".png");
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
    m_pimpl->viz.environment().lightViz("sun").setDirection(lightDirection);
    m_pimpl->textureForScreenshots->environment().lightViz("sun").setDirection(lightDirection);

    return true;

}

void Visualizer::visualizeWorldFrame(bool visualizeWorldFrame)
{
    m_pimpl->viz.environment().setElementVisibility("world_frame", visualizeWorldFrame);
    m_pimpl->textureForScreenshots->environment().setElementVisibility("world_frame", visualizeWorldFrame);
}

bool Visualizer::visualizeMPCStatesAndSaveAnimation(const std::vector<State> &states, std::function<iDynTree::Position (const State &)> stateCameraControl,
                                                    const std::vector<std::vector<State> > &fullStates, std::function<iDynTree::Position (const State &)> fullStateCameraControl,
                                                    const std::string &workingFolder, const std::string &fileName, const std::string &fileExtension, double speedUp, double endTime)
{
    if (!(m_pimpl->viz.getNrOfVisualizedModels())) {
        std::cerr << "[ERROR][Visualizer::visualizeMPCStatesAndSaveAnimation] First you have to load a model." << std::endl;
        return false;
    }

    if (states.size() == 0 || fullStates.size() == 0)
    {
        return true;
    }

    if (states.size() < fullStates.size())
    {
        std::cerr << "[ERROR][Visualizer::visualizeMPCStatesAndSaveAnimation] states must be equal or bigger than fullStates." << std::endl;
        return false;
    }

    unsigned int statesIncrement = (states.size() - 1) / fullStates.size(); //The first state is the initial one
    unsigned int numberOfFrames = fullStates.size() * (fullStates.front().size() + (statesIncrement - 1));

    unsigned int digits = static_cast<unsigned int>(std::floor(std::log10(numberOfFrames) + 1));

    size_t fullStateIndex = 0;
    size_t stateIndex = 1; //The first state is the initial one
    size_t frameIndex = 0;
    iDynTree::ILabel& speeduplabel = m_pimpl->viz.getLabel("speedup");
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2) << speedUp << "X";
    speeduplabel.setText(ss.str());
    speeduplabel.setPosition(iDynTree::Position(0 , states.front().comPosition[1] + 1.0, states.front().comPosition[2] + 1.0));

    DynamicalPlanner::State initialState = states.front();

    while (stateIndex < states.size() && fullStateIndex < fullStates.size() && (!(endTime < 0) || (states[stateIndex].time <= endTime))) {
        for (size_t i = 0; i < fullStates[fullStateIndex].size(); ++i)
        {
            speeduplabel.setVisible(false);
            m_pimpl->viz.camera().setPosition(stateCameraControl(initialState));
            m_pimpl->updateViz(initialState);
            m_pimpl->textureForScreenshots->setSubDrawArea(0, 0,
                                                           m_pimpl->textureForScreenshots->width()/2, m_pimpl->textureForScreenshots->height());
            m_pimpl->viz.subDraw(0, 0, m_pimpl->viz.width()/2, m_pimpl->viz.height());

            speeduplabel.setVisible(true);
            m_pimpl->viz.camera().setPosition(fullStateCameraControl(initialState));
            m_pimpl->updateViz(fullStates[fullStateIndex][i]);
            m_pimpl->textureForScreenshots->setSubDrawArea(m_pimpl->textureForScreenshots->width()/2, 0,
                                                           m_pimpl->textureForScreenshots->width()/2, m_pimpl->textureForScreenshots->height());
            m_pimpl->viz.subDraw(m_pimpl->viz.width()/2, 0, m_pimpl->viz.width()/2, m_pimpl->viz.height());

            m_pimpl->viz.draw();
            //Using a texture as the dimension should be always the same and not depending on the window size
            m_pimpl->textureForScreenshots->drawToFile(workingFolder + "/" + fileName + "_img_" + std::string(digits - std::to_string(frameIndex).size(), '0') + std::to_string(frameIndex) + ".png");
            frameIndex++;
        }

        for (size_t i = 0; i < statesIncrement; i++)
        {
            speeduplabel.setVisible(false);
            m_pimpl->viz.camera().setPosition(stateCameraControl(initialState));
            m_pimpl->updateViz(states[stateIndex]);
            m_pimpl->textureForScreenshots->setSubDrawArea(0, 0,
                                                           m_pimpl->textureForScreenshots->width()/2, m_pimpl->textureForScreenshots->height());
            m_pimpl->viz.subDraw(0, 0, m_pimpl->viz.width()/2, m_pimpl->viz.height());

            speeduplabel.setVisible(true);
            m_pimpl->viz.camera().setPosition(fullStateCameraControl(initialState));
            m_pimpl->updateViz(fullStates[fullStateIndex].back());
            m_pimpl->textureForScreenshots->setSubDrawArea(m_pimpl->textureForScreenshots->width()/2, 0,
                                                           m_pimpl->textureForScreenshots->width()/2, m_pimpl->textureForScreenshots->height());
            m_pimpl->viz.subDraw(m_pimpl->viz.width()/2, 0, m_pimpl->viz.width()/2, m_pimpl->viz.height());

            m_pimpl->viz.draw();
            //Using a texture as the dimension should be always the same and not depending on the window size
            m_pimpl->textureForScreenshots->drawToFile(workingFolder + "/" + fileName + "_img_" + std::string(digits - std::to_string(frameIndex).size(), '0') + std::to_string(frameIndex) + ".png");
            frameIndex++;

            initialState = states[stateIndex];
            stateIndex++;
        }
        fullStateIndex++;

        std::cout << "[INFO][Visualizer::visualizeState] Visualizing state " << stateIndex << "/" << states.size()
                  << ", iteration " << fullStateIndex << "/" << fullStates.size() << "." << std::endl;
    }

    speeduplabel.setVisible(false);

    m_pimpl->viz.camera().setPosition(m_pimpl->defaultCameraPosition);

    m_pimpl->viz.camera().setTarget(m_pimpl->defaultCameraTarget);

    std::string fps = std::to_string(static_cast<int>(std::round(stateIndex/states[stateIndex-1].time * speedUp))); //x2 speed

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

