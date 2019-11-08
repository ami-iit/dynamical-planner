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
#include <iDynTree/TimeRange.h>
#include <memory>
#include <string>
#include <vector>

namespace DynamicalPlanner {

    enum class ComplementarityType {
        HyperbolicSecantControl,
        Dynamical,
        Classical,
        HyperbolicTangentInDynamics
    };

    class Solver;
    class Settings;

    typedef struct {
        //Timings
        double minimumDt; //in seconds
        double maximumDt; //in seconds
        double controlPeriod; //in seconds
        double horizon; //in seconds
        double activeControlPercentage; //if 1, feet and forces can be changed for the full horizon. If 0 feet and forces are kept to the initial value

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


        //Complementarity Constraints
        ComplementarityType complementarity;
        //  - HyperbolicSecantControl
        iDynTree::Vector3 forceMaximumDerivative;
        double normalForceDissipationRatio;
        double normalForceHyperbolicSecantScaling;
        // - Dynamical Complementarity Constraint
        double complementarityDissipation;
        // - Classical Complementarity Constraint
        double classicalComplementarityTolerance;


        //ContactFrictionConstraint
        double frictionCoefficient;

        //ContactVelocityControlConstraints
        iDynTree::Vector3 velocityMaximumDerivative;
        double planarVelocityHyperbolicTangentScaling;
        double normalVelocityHyperbolicSecantScaling;
        bool contactVelocityControlConstraintsAsSeparateConstraints;



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
        std::vector<std::pair<double, double>> jointsVelocityLimits;
        bool constrainTargetCoMPosition;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> targetCoMPositionTolerance;
        iDynTree::optimalcontrol::TimeRange constrainTargetCoMPositionRange;



        //Costs
        //CoM cost
        bool comCostActive;
        double comCostOverallWeight;
        iDynTree::optimalcontrol::TimeRange comCostActiveRange;
        iDynTree::VectorDynSize comWeights;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMTrajectory;

        //CoM Velocity cost
        bool comVelocityCostActive;
        double comVelocityCostOverallWeight;
        iDynTree::VectorDynSize comVelocityWeights;
        iDynTree::optimalcontrol::TimeRange comVelocityCostActiveRange;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMVelocityTrajectory;

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

        //Swing height (each contact point has a different cost with same settings)
        bool swingCostActive;
        double swingCostOverallWeight;
        double desiredSwingHeight;
        iDynTree::Vector3 swingCostWeights;

        //Phantom forces aka forces when point is not in contact (each contact point has a different cost with same settings)
        bool phantomForcesCostActive;
        double phantomForcesCostOverallWeight;

        //Mean Position of the contact points
        bool meanPointPositionCostActive;
        double meanPointPositionCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> meanPointPositionCostTimeVaryingWeight;
        iDynTree::optimalcontrol::TimeRange meanPointPositionCostActiveRange;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> desiredMeanPointPosition;

        //Left foot yaw cost
        bool leftFootYawCostActive;
        double leftFootYawCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> desiredLeftFootYaw;

        //Right foot yaw cost
        bool rightFootYawCostActive;
        double rightFootYawCostOverallWeight;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> desiredRightFootYaw;

        //Right foot yaw cost
        bool feetDistanceCostActive;
        double feetDistanceCostOverallWeight;

        //Joint Velocity for postural cost regularization
        //(this cost share the same desired joints value and weights of joints regularization)
        bool jointsVelocityForPosturalCostActive;
        double jointsVelocityForPosturalCostOverallWeight;

        //Complementarity cost
        bool complementarityCostActive;
        double complementarityCostOverallWeight;

        //Base position cost
        bool basePositionCostActive;
        double basePositionCostOverallWeight;
        iDynTree::optimalcontrol::TimeRange basePositionCostActiveRange;
        iDynTree::VectorDynSize basePositionCostWeights;
        std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredBasePositionTrajectory;

        //Frame orientation (this uses the same frame of the frame orientation cost
        bool frameAngularVelocityCostActive;
        double frameAngularVelocityCostOverallWeight;
        double rotationalPIDgain;

        //Regularizations
        bool useConstraintsHessianRegularization;
        double constraintsHessianRegularization;
        bool useCostsHessianRegularization;
        double costsHessianRegularization;

    } SettingsStruct;
}

class DynamicalPlanner::Settings {

    DynamicalPlanner::SettingsStruct m_settings;
    bool m_isValid;

public:
    Settings();

    Settings(DynamicalPlanner::SettingsStruct& inputSettings);

    bool setFromStruct(const SettingsStruct &inputSettings);

    bool isValid() const;

    const SettingsStruct& getSettings() const;

    static DynamicalPlanner::SettingsStruct Defaults(const iDynTree::Model& newModel);

};

#endif // DPLANNER_SETTINGS_H
