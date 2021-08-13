/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Logger.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <matioCpp/matioCpp.h>
#include <cassert>

using namespace DynamicalPlanner;

template <typename Vector>
inline void addToLogger(std::unordered_map<std::string, matioCpp::MultiDimensionalArray<double>>& loggedVariables, const std::string& variableName, const Vector& value, size_t column) {

    matioCpp::MultiDimensionalArray<double>& var = loggedVariables[variableName];

    Eigen::Map<Eigen::MatrixXd> map(var.data(), var.dimensions()[0], var.dimensions()[1]);
    map.col(column) = iDynTree::toEigen(value);
}

void addToLogger(std::unordered_map<std::string, matioCpp::MultiDimensionalArray<double>>& loggedVariables, const std::string& variableName, double value, size_t column) {

    loggedVariables[variableName]({0, column}) = value;
}

template <typename Type>
inline void setMatioStructField(matioCpp::Struct& structToFill, const std::string& name, const Type& scalarValue)
{
    structToFill.setField(matioCpp::Element<Type>(name, scalarValue));
}

template <>
inline void setMatioStructField<bool>(matioCpp::Struct& structToFill, const std::string& name, const bool& scalarValue)
{
    structToFill.setField(matioCpp::Element<matioCpp::Logical>(name, scalarValue));
}

template <>
inline void setMatioStructField<std::string>(matioCpp::Struct& structToFill, const std::string& name, const std::string& scalarValue)
{
    structToFill.setField(matioCpp::String(name, scalarValue));
}

template <>
inline void setMatioStructField<iDynTree::Vector3>(matioCpp::Struct& structToFill, const std::string& name, const iDynTree::Vector3& scalarValue)
{
    structToFill.setField(matioCpp::Vector<double>(name, scalarValue));
}

template <>
inline void setMatioStructField<iDynTree::VectorDynSize>(matioCpp::Struct& structToFill, const std::string& name, const iDynTree::VectorDynSize& scalarValue)
{
    structToFill.setField(matioCpp::Vector<double>(name, scalarValue));
}

template <>
inline void setMatioStructField<std::vector<double>>(matioCpp::Struct& structToFill, const std::string& name, const std::vector<double>& scalarValue)
{
    structToFill.setField(matioCpp::Vector<double>(name, scalarValue));
}

#define SET_MATIO_STRUCT_FIELD(structToFill, inputStruct, field) setMatioStructField(structToFill, #field, inputStruct.field)

matioCpp::Struct populateSettingsStruct(const SettingsStruct& settings)
{
    matioCpp::Struct settingsVar("settings");

    SET_MATIO_STRUCT_FIELD(settingsVar, settings, minimumDt);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, maximumDt);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, controlPeriod);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, horizon);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, activeControlPercentage);
    settingsVar.setField(matioCpp::String("robotModel", settings.robotModel.toString()));
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, gravity);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, updateTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, floatingBaseName);

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

    SET_MATIO_STRUCT_FIELD(settingsVar, settings, leftFrameName);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, rightFrameName);

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

    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceMaximumDerivative);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, normalForceDissipationRatio);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, normalForceHyperbolicSecantScaling);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, complementarityDissipation);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, dynamicComplementarityUpperBound);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, classicalComplementarityTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, frictionCoefficient);

    switch (settings.planarComplementarity)
    {
    case PlanarComplementarityType::Classical:
        settingsVar.setField(matioCpp::String("planarComplementarity", "Classical"));
        break;
    case PlanarComplementarityType::HyperbolicTangentInDynamics:
        settingsVar.setField(matioCpp::String("planarComplementarity", "HyperbolicTangentInDynamics"));
        break;
    case PlanarComplementarityType::HyperbolicTangentInequality:
        settingsVar.setField(matioCpp::String("planarComplementarity", "HyperbolicTangentInequality"));
        break;
    }

    SET_MATIO_STRUCT_FIELD(settingsVar, settings, velocityMaximumDerivative);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, planarVelocityHyperbolicTangentScaling);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, normalVelocityHyperbolicSecantScaling);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, classicalPlanarComplementarityTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, indexOfLateralDirection);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, minimumFeetDistance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, referenceFrameNameForFeetDistance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, otherFrameNameForFeetDistance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, feetMaximumRelativeHeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comPositionConstraintTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, centroidalMomentumConstraintTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, quaternionModulusConstraintTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, pointPositionConstraintTolerance);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, minimumCoMHeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, maximumAngularMomentum);

    std::vector<matioCpp::Variable> bounds;
    for (auto& bound : settings.jointsLimits)
    {
        matioCpp::Struct boundStruct("bound");
        boundStruct.setField(matioCpp::Element<double>("min", bound.first));
        boundStruct.setField(matioCpp::Element<double>("max", bound.second));
        bounds.push_back(boundStruct);
    }
    settingsVar.setField(matioCpp::CellArray("jointsLimits", {bounds.size(), 1}, bounds));

    bounds.clear();
    for (auto& bound : settings.jointsVelocityLimits)
    {
        matioCpp::Struct boundStruct("bound");
        boundStruct.setField(matioCpp::Element<double>("min", bound.first));
        boundStruct.setField(matioCpp::Element<double>("max", bound.second));
        bounds.push_back(boundStruct);
    }
    settingsVar.setField(matioCpp::CellArray("jointsVelocityLimits", {bounds.size(), 1}, bounds));

    SET_MATIO_STRUCT_FIELD(settingsVar, settings, constrainTargetCoMPosition);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comVelocityCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comVelocityCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, comVelocityWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, frameCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, frameCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, frameForOrientationCost);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceMeanCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceMeanCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceRatioCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceRatioCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, desiredLeftRatios);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, desiredRightRatios);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsRegularizationCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsRegularizationCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsRegularizationWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsVelocityCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsVelocityCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsVelocityCostWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, staticTorquesCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, staticTorquesCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, staticTorquesCostWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceDerivativeCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceDerivativesCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, forceDerivativeWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, pointAccelerationCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, pointAccelerationCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, pointAccelerationWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, swingCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, swingCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, desiredSwingHeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, swingCostWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, phantomForcesCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, phantomForcesCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, meanPointPositionCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, meanPointPositionCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, leftFootYawCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, leftFootYawCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, rightFootYawCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, rightFootYawCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, feetDistanceCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, feetDistanceCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsVelocityForPosturalCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, jointsVelocityForPosturalCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, complementarityCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, complementarityCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, basePositionCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, basePositionCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, basePositionCostWeights);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, baseQuaternionCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, baseQuaternionCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, frameAngularVelocityCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, frameAngularVelocityCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, rotationalPIDgain);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, baseQuaternionVelocityCostActive);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, baseQuaternionVelocityCostOverallWeight);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, useConstraintsHessianRegularization);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, constraintsHessianRegularization);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, useCostsHessianRegularization);
    SET_MATIO_STRUCT_FIELD(settingsVar, settings, costsHessianRegularization);

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

