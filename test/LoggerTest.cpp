/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Logger.h>
#include <iDynTree/ModelIO/ModelLoader.h>
#include <URDFdir.h>
#include <FolderPath.h>
#include <chrono>

int main()
{
    std::vector<DynamicalPlanner::State> states(4, DynamicalPlanner::State(23, 4));
    std::vector<DynamicalPlanner::Control> controls(5, DynamicalPlanner::Control(23, 4));

    iDynTree::ModelLoader modelLoader;
    modelLoader.loadModelFromFile(getAbsModelPath("iCubGenova04.urdf"));
    DynamicalPlanner::SettingsStruct settings = DynamicalPlanner::Settings::Defaults(modelLoader.model());

    states[0].time = 1.0;
    controls[0].time = 2.0;

    auto timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm timeStruct;
    timeStruct = *std::localtime(&timeNow);
    std::ostringstream timeString;

    timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    timeStruct = *std::localtime(&timeNow);
    timeString << timeStruct.tm_year + 1900 << "-" << timeStruct.tm_mon << "-";
    timeString << timeStruct.tm_mday << "_" << timeStruct.tm_hour << "_" << timeStruct.tm_min;
    timeString << "_" << timeStruct.tm_sec;

    DynamicalPlanner::Logger::saveSolutionVectorsToFile(getAbsDirPath("SavedVideos") + "/log-" + timeString.str() + ".mat" , settings, states, controls);

    return EXIT_SUCCESS;
}
