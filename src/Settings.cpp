/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Settings.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iostream>
#include <sstream>

using namespace DynamicalPlanner;


void printError(const std::string& message) {
    std::cerr << "[ERROR][Settings::setFromStruct] "<< message <<std::endl;
}

int checkError(bool error, const std::string& errorMessage) {
    if (error) {
        printError(errorMessage);
    }
    return error;
}


Settings::Settings()
    : m_isValid(false)
{ }

Settings::Settings(SettingsStruct &inputSettings)
{
    m_isValid = setFromStruct(inputSettings);
}

bool Settings::setFromStruct(const SettingsStruct &inputSettings)
{
    std::ostringstream errorlabel;
    errorlabel << "[ERROR][Settings::setFromStruct] ";

    int errors = 0;

    errors += checkError(inputSettings.minimumDt < 0, "The minimumDt is supposed to be non-negative.");
    errors += checkError(inputSettings.maximumDt < 0, "The maximumDt is supposed to be non-negative.");
    errors += checkError(inputSettings.maximumDt <= 2*inputSettings.minimumDt,
                         "The minimumDt is supposed to be lower than half of the maximumDt.");
    errors += checkError(inputSettings.controlPeriod < 0, "The controlPeriod is supposed to be non-negative.");
    errors += checkError(inputSettings.horizon < 0, "The horizon is supposed to be non-negative.");

    //errors += checkError(inputSettings.updateTolerance < 0, "The updateTolerance is supposed to be non-negative.");
    errors += checkError(!(inputSettings.robotModel.isLinkNameUsed(inputSettings.floatingBaseName)),
                         "The floating base name has to refer to a link in the model.");
    errors += checkError(inputSettings.leftPointsPosition.size() != inputSettings.rightPointsPosition.size(),
                         "It is expected to have the same number of points on the left and right foot.");
    errors += checkError(!(inputSettings.robotModel.isFrameNameUsed(inputSettings.leftFrameName)),
                         "The leftFrameName does not appear in the model.");
    errors += checkError(!(inputSettings.robotModel.isFrameNameUsed(inputSettings.rightFrameName)),
                         "The rightFrameName does not appear in the model.");
    errors += checkError(inputSettings.complementarityDissipation < 0, "The complementarityDissipation is expected to be non-negative.");
    errors += checkError(inputSettings.frictionCoefficient < 0, "The frictionCoefficient is expected to be non-negative.");
    errors += checkError(inputSettings.indexOfLateralDirection > 2, "The indexOfLateralDirection is expected to be in the range [0, 2].");
    errors += checkError(inputSettings.minimumFeetDistance < 0, "The minimumDistance is expected to be non-negative.");
    errors += checkError(!(inputSettings.robotModel.isFrameNameUsed(inputSettings.referenceFrameNameForFeetDistance)),
                         "The referenceFrameName does not appear in the model.");
    errors += checkError(!(inputSettings.robotModel.isFrameNameUsed(inputSettings.otherFrameNameForFeetDistance)),
                         "The otherFrameName does not appear in the model.");
    errors += checkError(inputSettings.comPositionConstraintTolerance < 0, "The comPositionConstraintTolerance is expected to be non-negative.");
    errors += checkError(inputSettings.centroidalMomentumConstraintTolerance < 0,
                         "The centroidalMomentumConstraintTolerance is expected to be non-negative.");
    errors += checkError(inputSettings.quaternionModulusConstraintTolerance < 0,
                         "The quaternionModulusConstraintTolerance is expected to be non-negative.");
    errors += checkError(inputSettings.pointPositionConstraintTolerance < 0,
                         "The pointPositionConstraintTolerance is expected to be non-negative.");
    errors += checkError(inputSettings.minimumCoMHeight < 0, "The minimumCoMHeight is expected to be non-negative.");
    errors += checkError(inputSettings.jointsLimits.size() != inputSettings.robotModel.getNrOfDOFs(),
                         "The jointsLimits size differs from the number of joints.");
    errors += checkError(inputSettings.jointsVelocityLimits.size() != inputSettings.robotModel.getNrOfDOFs(),
                         "The jointsVelocityLimits size differs from the number of joints.");
    for (size_t j = 0; j < inputSettings.jointsLimits.size(); ++j) {
        errors += checkError(inputSettings.jointsLimits[j].first > inputSettings.jointsLimits[j].second,
                             "The lower limit of joint " + std::to_string(j) + " is bigger than its upper limit.");

        errors += checkError(inputSettings.jointsVelocityLimits[j].first > inputSettings.jointsVelocityLimits[j].second,
                             "The lower velocity limit of joint " + std::to_string(j) + " is bigger than its upper limit.");
    }

    if (inputSettings.comCostActive) {
        errors += checkError(inputSettings.comWeights.size() != 3, "The comWeight vector is expected to be three dimensional.");
        errors += checkError(!inputSettings.desiredCoMTrajectory, "The desiredCoMTrajectory is empty.");
    }

    if (inputSettings.frameCostActive) {
        errors += checkError(!(inputSettings.robotModel.isFrameNameUsed(inputSettings.frameForOrientationCost)),
                             "The frameForOrientationCost does not appear in the model.");
        errors += checkError(!inputSettings.desiredRotationTrajectory, "The desiredRotationTrajectory is empty.");
    }

    if (inputSettings.jointsRegularizationCostActive) {
        errors += checkError(!(inputSettings.desiredJointsTrajectory), "The desiredJointsTrajectory is empty.");
        errors += checkError(inputSettings.jointsRegularizationWeights.size() != inputSettings.robotModel.getNrOfDOFs(),
                             "The jointRegularizationWeights size does not match the robot DoFs.");
    }

    if (inputSettings.staticTorquesCostActive) {
        errors += checkError(inputSettings.staticTorquesCostWeights.size() != inputSettings.robotModel.getNrOfDOFs(),
                             "The staticTorquesCostWeights size does not match the robot DoFs.");
    }

    if (inputSettings.forceDerivativeCostActive) {
        errors += checkError(!(inputSettings.desiredRotationTrajectory), " The desiredRotationTrajectory is empty.");
        errors += checkError(inputSettings.forceDerivativeWeights.size() != 3, "The forceDerivativeWeights is expected to be three dimensional.");
    }

    if (inputSettings.pointAccelerationCostActive) {
        errors += checkError(!(inputSettings.desiredPointAccelerationTrajectory), " The desiredPointVelocityTrajectory is empty.");
        errors += checkError(inputSettings.pointAccelerationWeights.size() != 3, "The pointVelocityWeights is expected to be three dimensional.");
    }

    checkError(errors > 0, "The were errors when importing the settings struct. The settings will not be updated.");

    if (errors == 0) {
        m_settings = inputSettings;
        m_isValid = true;
    }

    return (errors == 0);
}

bool Settings::isValid() const
{
    return m_isValid;
}

const SettingsStruct &Settings::getSettings() const
{
    return m_settings;
}

SettingsStruct Settings::Defaults(const iDynTree::Model &newModel)
{
    SettingsStruct defaults;

    defaults.minimumDt = 0.01;
    defaults.maximumDt = 0.1;
    defaults.controlPeriod = 0.01;
    defaults.horizon = 1.0;

        // SharedKinDyn
    defaults.robotModel = newModel;
    defaults.gravity.zero();
    defaults.gravity(2) = -9.81;

    defaults.updateTolerance = 1e-10;
    defaults.floatingBaseName = "root_link";

    //Contact points infos
    defaults.leftPointsPosition.resize(4);
    iDynTree::toEigen(defaults.leftPointsPosition[0]) << 0.125, -0.04, 0.0;
    iDynTree::toEigen(defaults.leftPointsPosition[1]) << 0.125,  0.04, 0.0;
    iDynTree::toEigen(defaults.leftPointsPosition[2]) << -0.063,  0.04, 0.0;
    iDynTree::toEigen(defaults.leftPointsPosition[3]) << 0.063, -0.04, 0.0;
    defaults.rightPointsPosition.resize(4);
    iDynTree::toEigen(defaults.rightPointsPosition[0]) << 0.125,  0.04, 0.0;
    iDynTree::toEigen(defaults.rightPointsPosition[1]) << 0.125, -0.04, 0.0;
    iDynTree::toEigen(defaults.rightPointsPosition[2]) << -0.063, -0.04, 0.0;
    iDynTree::toEigen(defaults.rightPointsPosition[3]) << 0.063,  0.04, 0.0;
    defaults.leftFrameName = "l_sole";
    defaults.rightFrameName = "r_sole";

    //ContactForceControlConstraints
    iDynTree::toEigen(defaults.forceMaximumDerivative).setConstant(10.0);
    defaults.normalForceDissipationRatio = 10.0;
    defaults.normalForceHyperbolicSecantScaling = 300.0;

    //Dynamical Complementarity Constraint
    defaults.complementarityDissipation = 10.0;

    //ContactFrictionConstraint
    defaults.frictionCoefficient = 0.3;

    //ContactVelocityControlConstraints
    iDynTree::toEigen(defaults.velocityMaximumDerivative).setConstant(10.0);
    defaults.planarVelocityHyperbolicTangentScaling = 50.0;
    defaults.normalVelocityHyperbolicSecantScaling = 1.0;

    //Feet lateral distance constraint
    defaults.indexOfLateralDirection = 1;
    defaults.minimumFeetDistance = 0.1;
    defaults.referenceFrameNameForFeetDistance = "r_sole";
    defaults.otherFrameNameForFeetDistance = "l_sole";

    //Equality constraints tolerances
    defaults.comPositionConstraintTolerance = 1e-3;
    defaults.centroidalMomentumConstraintTolerance = 1e-2;
    defaults.quaternionModulusConstraintTolerance = 1e-4;
    defaults.pointPositionConstraintTolerance = 1e-3;

    //Bounds
    defaults.minimumCoMHeight = 0.2;

    unsigned int n = static_cast<unsigned int>(newModel.getNrOfDOFs());

    std::pair<double, double> singleJointLimit;

    //for each joint, ask the limits
    for (iDynTree::JointIndex jointIdx = 0; jointIdx < static_cast<int>(n); ++jointIdx) {
        iDynTree::IJointConstPtr joint = newModel.getJoint(jointIdx);
        //if the joint does not have limits skip it
        if (!joint->hasPosLimits())
            continue;
        //for each DoF modelled by the joint get the limits
        for (unsigned dof = 0; dof < joint->getNrOfDOFs(); ++dof) {
            if (!joint->getPosLimits(dof, singleJointLimit.first, singleJointLimit.second))
                continue;

            defaults.jointsLimits.push_back(singleJointLimit);
        }
    }

    std::pair<double, double> defaultJointsVelocityLimits(-1E19, 1E19);
    defaults.jointsVelocityLimits.resize(newModel.getNrOfDOFs(), defaultJointsVelocityLimits);

    //Costs
    //CoM cost
    defaults.comCostActive = true;
    defaults.comCostOverallWeight = 1.0;
    defaults.comWeights.resize(3);
    iDynTree::toEigen(defaults.comWeights).setConstant(1.0);
    iDynTree::VectorDynSize desiredCoM(3);
    desiredCoM.zero();
    desiredCoM(2) = 0.5;
    defaults.desiredCoMTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredCoM);

    //Frame orientation
    defaults.frameCostActive = true;
    defaults.frameCostOverallWeight = 0.5;
    defaults.frameForOrientationCost = "neck_2";
    iDynTree::Rotation desiredRotation(0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    defaults.desiredRotationTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantRotation>(desiredRotation);

    //Force mean (each contact point has a different cost with same settings)
    defaults.forceMeanCostActive = true;
    defaults.forceMeanCostOverallWeight = 0.01;


    //Joints regularization
    defaults.jointsRegularizationCostActive = true;
    defaults.jointsRegularizationCostOverallWeight = 0.001;
    iDynTree::VectorDynSize desiredJointsConfiguration(n);
    desiredJointsConfiguration.zero();
    defaults.desiredJointsTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredJointsConfiguration);
    defaults.jointsRegularizationWeights.resize(n);
    iDynTree::toEigen(defaults.jointsRegularizationWeights).setConstant(1.0);

    //Joints velocity regularization
    defaults.jointsVelocityCostActive = true;
    defaults.jointsVelocityCostOverallWeight = 0.0001;
    defaults.desiredJointsVelocityTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredJointsConfiguration);;
    defaults.jointsVelocityCostWeights = defaults.jointsRegularizationWeights;

    //Static torques
    defaults.staticTorquesCostActive = true;
    defaults.staticTorquesCostOverallWeight = 0.01;
    defaults.staticTorquesCostWeights = defaults.jointsRegularizationWeights;

    //Force derivative cost (each contact point has a different cost with same settings)
    defaults.forceDerivativeCostActive = true;
    defaults.forceDerivativesCostOverallWeight = 0.0001;
    iDynTree::VectorDynSize desiredForceDerivative(3);
    desiredForceDerivative.zero();
    defaults.desiredForceDerivativeTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredForceDerivative);
    defaults.forceDerivativeWeights.resize(3);
    iDynTree::toEigen(defaults.forceDerivativeWeights).setConstant(1.0);

    //Point velocity cost (each contact point has a different cost with same settings)
    defaults.pointAccelerationCostActive = true;
    defaults.pointAccelerationCostOverallWeight = 0.0001;
    iDynTree::VectorDynSize desiredPointVelocity(3);
    desiredPointVelocity.zero();
    defaults.desiredPointAccelerationTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(desiredPointVelocity);
    defaults.pointAccelerationWeights.resize(3);
    iDynTree::toEigen(defaults.pointAccelerationWeights).setConstant(1.0);

    return defaults;
}
