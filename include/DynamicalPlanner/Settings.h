/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SETTINGS_H
#define DPLANNER_SETTINGS_H

#include <iDynTree/Model/Model.h>
#include <iDynTree/Core/Position.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/TimeVaryingObject.h>
#include <memory>
#include <string>
#include <vector>

namespace DynamicalPlanner {
    class Solver;
    class Settings;

    typedef struct {
        //Timings
        double minimumDt; //in seconds
        double maximumDt; //in seconds
        double controlPeriod; //in seconds
        double horizon; //in seconds

        // SharedKinDyn
        iDynTree::Model robotModel;
        iDynTree::Vector3 gravity;
        double updateTolerance;
        std::string floatingBaseName;

        //Contact points infos
        std::vector<iDynTree::Position> leftPointsPosition;
        std::vector<iDynTree::Position> rightPointsPosition;
        std::string leftFrameName;
        std::string rightFrameName;

        //ContactForceControlConstraints
        iDynTree::Vector3 forceMaximumDerivative;
        iDynTree::Vector3 forceDissipationRatio;
        double forceHyperbolicSecantScaling;

        //ContactFrictionConstraint
        double frictionCoefficient;

        //ContactVelocityControlConstraints
        iDynTree::Vector3 velocityMaximumDerivative;
        iDynTree::Vector3 velocityDissipationRatio;
        double velocityHyperbolicSecantScalingXY; //scales the position along z
        double velocityHyperbolicSecantScalingZ; //scales the force along z

        //Feet lateral distance constraint
        unsigned int indexOfLateralDirection;
        double minimumFeetDistance;
        std::string referenceFrameNameForFeetDistance;
        std::string otherFrameNameForFeetDistance;

        //Equality constraints tolerances
        double comPositionConstraintTolerance;
        double centroidalMomentumConstraintTolerance;
        double quaternionModulusConstraintTolerance;
        double pointPositionConstraintTolerance;

        //Bounds
        double minimumCoMHeight;
        std::vector<std::pair<double, double>> jointsLimits;

        //Costs
        //CoM cost
        bool comCostActive;
        double comCostOverallWeight;
        iDynTree::VectorDynSize comWeights;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMTrajectory;

        //Frame orientation
        bool frameCostActive;
        double frameCostOverallWeight;
        std::string frameForOrientationCost;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingRotation> desiredRotationTrajectory;

        //Force mean (each contact point has a different cost with same settings)
        bool forceMeanCostActive;
        double forceMeanCostOverallWeight;

        //Joint regularization
        bool jointsRegularizationCostActive;
        double jointsRegularizationCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredJointsTrajectory;
        iDynTree::VectorDynSize jointsRegularizationWeights;

        //Joint regularization
        bool jointsVelocityCostActive;
        double jointsVelocityCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredJointsVelocityTrajectory;
        iDynTree::VectorDynSize jointsVelocityCostWeights;

        //Static torques
        bool staticTorquesCostActive;
        double staticTorquesCostOverallWeight;
        iDynTree::VectorDynSize staticTorquesCostWeights;

        //Force derivative cost (each contact point has a different cost with same settings)
        bool forceDerivativeCostActive;
        double forceDerivativesCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredForceDerivativeTrajectory;
        iDynTree::VectorDynSize forceDerivativeWeights;

        //Point acceleration cost (each contact point has a different cost with same settings)
        bool pointAccelerationCostActive;
        double pointAccelerationCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredPointAccelerationTrajectory;
        iDynTree::VectorDynSize pointAccelerationWeights;

    } SettingsStruct;
}

class DynamicalPlanner::Settings {

    DynamicalPlanner::SettingsStruct m_settings;
    bool m_isValid;

public:
    Settings();

    Settings(DynamicalPlanner::SettingsStruct& inputSettings);

    bool setFromStruct(DynamicalPlanner::SettingsStruct& inputSettings);

    bool isValid() const;

    const SettingsStruct& getSettings() const;

    static DynamicalPlanner::SettingsStruct Defaults(const iDynTree::Model& newModel);

};

#endif // DPLANNER_SETTINGS_H
