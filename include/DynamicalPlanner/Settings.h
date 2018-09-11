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
        double minimumDistance;
        std::string referenceFrameName;
        std::string otherFrameName;

        //Equality constraints tolerances
        double comPositionConstraintTolerance;
        double centroidalMomentumConstraintTolerance;
        double quaternionModulusConstraintTolerance;
        double pointPositionConstraintTolerance;

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
        double jjointsVelocityCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredJointsVelocityTrajectory;
        iDynTree::VectorDynSize jointsVelocityCostWeights;

        //Static torques
        bool staticTorquesCostActive;
        double staticTorquesCostOverallWeight;
        iDynTree::VectorDynSize staticTorquesCostWeights;

        //Force derivative cost (each contact point has a different cost with same settings)
        bool forceDerivativeCostActive;
        bool forceDerivativesCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredForceDerivativeTrajectory;
        iDynTree::VectorDynSize forceDerivativeWeights;

        //Point velocity cost (each contact point has a different cost with same settings)
        bool pointVelocityCostActive;
        bool pointVelocityCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredPointVelocityTrajectory;
        iDynTree::VectorDynSize pointVelocityWeights;

    } SettingsStruct;
}

class DynamicalPlanner::Settings {

    friend class DynamicalPlanner::Solver;

    DynamicalPlanner::SettingsStruct m_settings;

public:
    Settings();

    bool setFromStruct(DynamicalPlanner::SettingsStruct& inputSettings);

    static DynamicalPlanner::SettingsStruct Defaults();

};

#endif // DPLANNER_SETTINGS_H
