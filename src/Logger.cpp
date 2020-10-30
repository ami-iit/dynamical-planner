/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Logger.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <matioCpp/matioCpp.h>

using namespace DynamicalPlanner;

template <typename Vector>
void addToLogger(std::unordered_map<std::string, matioCpp::MultiDimensionalArray<double>>& loggedVariables, const std::string& variableName, const Vector& value, size_t column) {

    matioCpp::MultiDimensionalArray<double>& var = loggedVariables[variableName];

    Eigen::Map<Eigen::MatrixXd> map(var.data(), var.dimensions()[0], var.dimensions()[1]);
    map.col(column) = iDynTree::toEigen(value);
}

void addToLogger(std::unordered_map<std::string, matioCpp::MultiDimensionalArray<double>>& loggedVariables, const std::string& variableName, double value, size_t column) {

    loggedVariables[variableName]({0, column}) = value;
}

matioCpp::Struct populateSettingsStruct(const SettingsStruct& settings)
{
    matioCpp::Struct settingsVar("settings");

    settingsVar.setField(matioCpp::Element<double>("minimumDt", settings.minimumDt));
    settingsVar.setField(matioCpp::Element<double>("maximumDt", settings.maximumDt));
    settingsVar.setField(matioCpp::Element<double>("controlPeriod", settings.controlPeriod));
    settingsVar.setField(matioCpp::Element<double>("horizon", settings.horizon));
    settingsVar.setField(matioCpp::Element<double>("activeControlPercentage", settings.activeControlPercentage));

    settingsVar.setField(matioCpp::String("robotModel", settings.robotModel.toString()));
    settingsVar.setField(matioCpp::Vector<double>("gravity", settings.gravity));
    settingsVar.setField(matioCpp::Element<double>("updateTolerance", settings.updateTolerance));
    settingsVar.setField(matioCpp::String("floatingBaseName", settings.floatingBaseName));

    std::vector<matioCpp::Variable> pointsPosition;
    for (auto& pos : settings.leftPointsPosition)
    {
        pointsPosition.push_back(matioCpp::Vector<double>("contactPoint", pos));
    }
    settingsVar.setField(matioCpp::CellArray("leftPointsPosition", {settings.leftPointsPosition.size(), 1}, pointsPosition));

    pointsPosition.clear();
    for (auto& pos : settings.rightPointsPosition)
    {
        pointsPosition.push_back(matioCpp::Vector<double>("contactPoint", pos));
    }
    settingsVar.setField(matioCpp::CellArray("rightPointsPosition", {settings.rightPointsPosition.size(), 1}, pointsPosition));
    settingsVar.setField(matioCpp::String("leftFrameName", settings.leftFrameName));
    settingsVar.setField(matioCpp::String("rightFrameName", settings.rightFrameName));

    switch (settings.complementarity)
    {
    case ComplementarityType::Classical:
        settingsVar.setField(matioCpp::String("complementarity", "Classical"));
        break;
    case ComplementarityType::Dynamical:
        settingsVar.setField(matioCpp::String("complementarity", "Dynamical"));
        break;
    case ComplementarityType::HyperbolicSecantInDynamics:
        settingsVar.setField(matioCpp::String("complementarity", "HyperbolicSecantInDynamics"));
        break;

    case ComplementarityType::HyperbolicSecantInequality:
        settingsVar.setField(matioCpp::String("complementarity", "HyperbolicSecantInequality"));
        break;
    }

    settingsVar.setField(matioCpp::Vector<double>("forceMaximumDerivative", settings.forceMaximumDerivative));
    settingsVar.setField(matioCpp::Element<double>("normalForceDissipationRatio", settings.normalForceDissipationRatio));
    settingsVar.setField(matioCpp::Element<double>("normalForceHyperbolicSecantScaling", settings.normalForceHyperbolicSecantScaling));
    settingsVar.setField(matioCpp::Element<double>("complementarityDissipation", settings.complementarityDissipation));
    settingsVar.setField(matioCpp::Element<double>("dynamicComplementarityUpperBound", settings.dynamicComplementarityUpperBound));
    settingsVar.setField(matioCpp::Element<double>("classicalComplementarityTolerance", settings.classicalComplementarityTolerance));


    return settingsVar;
}

void Logger::saveSolutionVectorsToFile(const std::string &matFileName, const SettingsStruct& settings, const std::vector<State> &states, const std::vector<Control> &controls, const std::vector<double> &computationalTime)
{
    std::unordered_map<std::string, matioCpp::MultiDimensionalArray<double>> loggedVariables;

    if (states.size()) {
        for (size_t point = 0; point < states.front().leftContactPointsState.size(); ++point) {
            std::string name = "leftPoint" + std::to_string(point) + "Force";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, states.size()});
            name = "leftPoint" + std::to_string(point) + "Position";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, states.size()});
        }

        for (size_t point = 0; point < states.front().rightContactPointsState.size(); ++point) {
            std::string name = "rightPoint" + std::to_string(point) + "Force";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, states.size()});
            name = "rightPoint" + std::to_string(point) + "Position";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, states.size()});
        }

        loggedVariables["momentumInCoM"] = matioCpp::MultiDimensionalArray<double>("momentumInCoM", {6, states.size()});
        loggedVariables["comPosition"] = matioCpp::MultiDimensionalArray<double>("comPosition", {3, states.size()});
        loggedVariables["basePosition"] = matioCpp::MultiDimensionalArray<double>("basePosition", {3, states.size()});
        loggedVariables["baseQuaternion"] = matioCpp::MultiDimensionalArray<double>("baseQuaternion", {4, states.size()});
        loggedVariables["jointsConfiguraion"] = matioCpp::MultiDimensionalArray<double>("jointsConfiguraion", {states.front().jointsConfiguration.size(), states.size()});
        loggedVariables["stateTime"] = matioCpp::MultiDimensionalArray<double>("stateTime", {1, states.size()});
    }

    if (controls.size()) {
        for (size_t point = 0; point < controls.front().leftContactPointsControl.size(); ++point) {
            std::string name = "leftPoint" + std::to_string(point) + "ForceControl";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, controls.size()});
            name = "leftPoint" + std::to_string(point) + "VelocityControl";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, controls.size()});
        }

        for (size_t point = 0; point < controls.front().rightContactPointsControl.size(); ++point) {
            std::string name = "rightPoint" + std::to_string(point) + "ForceControl";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, controls.size()});
            name = "rightPoint" + std::to_string(point) + "VelocityControl";
            loggedVariables[name] = matioCpp::MultiDimensionalArray<double>(name, {3, controls.size()});
        }

        loggedVariables["baseLinearVelocity"] = matioCpp::MultiDimensionalArray<double>("baseLinearVelocity", {3, controls.size()});
        loggedVariables["baseQuaternionDerivative"] = matioCpp::MultiDimensionalArray<double>("baseQuaternionDerivative", {4, controls.size()});
        loggedVariables["jointsVelocity"] = matioCpp::MultiDimensionalArray<double>("jointsVelocity", {controls.front().jointsVelocity.size(), controls.size()});
        loggedVariables["controlTime"] = matioCpp::MultiDimensionalArray<double>("controlTime", {1, controls.size()});

    }

    if (computationalTime.size())
    {
        loggedVariables["computationalTime"] = matioCpp::MultiDimensionalArray<double>("computationalTime", {1, computationalTime.size()});
    }

    for (size_t i = 0; i < states.size(); ++i) {
        for (size_t point = 0; point < states.front().leftContactPointsState.size(); ++point) {
            addToLogger(loggedVariables, "leftPoint" + std::to_string(point) + "Force", states[i].leftContactPointsState[point].pointForce, i);
            addToLogger(loggedVariables,"leftPoint" + std::to_string(point) + "Position", states[i].leftContactPointsState[point].pointPosition, i);
        }

        for (size_t point = 0; point < states.front().rightContactPointsState.size(); ++point) {
            addToLogger(loggedVariables, "rightPoint" + std::to_string(point) + "Force", states[i].rightContactPointsState[point].pointForce, i);
            addToLogger(loggedVariables,"rightPoint" + std::to_string(point) + "Position", states[i].rightContactPointsState[point].pointPosition, i);
        }

        addToLogger(loggedVariables, "momentumInCoM", states[i].momentumInCoM, i);
        addToLogger(loggedVariables, "comPosition", states[i].comPosition, i);
        addToLogger(loggedVariables, "basePosition", states[i].worldToBaseTransform.getPosition(), i);
        addToLogger(loggedVariables, "baseQuaternion", states[i].worldToBaseTransform.getRotation().asQuaternion(), i);
        addToLogger(loggedVariables, "jointsConfiguraion", states[i].jointsConfiguration, i);
        addToLogger(loggedVariables, "stateTime", states[i].time, i);

    }

    for (size_t i = 0; i < controls.size(); ++i) {
        for (size_t point = 0; point < controls.front().leftContactPointsControl.size(); ++point) {
            addToLogger(loggedVariables, "leftPoint" + std::to_string(point) + "ForceControl", controls[i].leftContactPointsControl[point].pointForceControl, i);
            addToLogger(loggedVariables, "leftPoint" + std::to_string(point) + "VelocityControl",
                        controls[i].leftContactPointsControl[point].pointVelocityControl, i);
        }

        for (size_t point = 0; point < controls.front().rightContactPointsControl.size(); ++point) {
            addToLogger(loggedVariables, "rightPoint" + std::to_string(point) + "ForceControl", controls[i].rightContactPointsControl[point].pointForceControl, i);
            addToLogger(loggedVariables, "rightPoint" + std::to_string(point) + "VelocityControl",
                        controls[i].rightContactPointsControl[point].pointVelocityControl, i);
        }

        addToLogger(loggedVariables, "baseLinearVelocity", controls[i].baseLinearVelocity, i);
        addToLogger(loggedVariables, "baseQuaternionDerivative", controls[i].baseQuaternionDerivative, i);
        addToLogger(loggedVariables, "jointsVelocity", controls[i].jointsVelocity, i);
        addToLogger(loggedVariables, "controlTime", controls[i].time, i);

    }

    for (size_t i = 0; i < computationalTime.size(); ++i)
    {
        addToLogger(loggedVariables, "computationalTime", computationalTime[i], i);
    }

    matioCpp::File file = matioCpp::File::Create(matFileName);

    file.write(loggedVariables.begin(), loggedVariables.end());
    file.write(populateSettingsStruct(settings));

}

