/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SOLVER_H
#define DPLANNER_SOLVER_H

#include <iDynTree/Optimizer.h>
#include <DynamicalPlanner/Settings.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/Transform.h>
#include <memory>

namespace DynamicalPlanner {
    class Solver;

    typedef struct {
        iDynTree::Vector3 pointForce;
        iDynTree::Vector3 pointVelocity;
        iDynTree::Vector3 pointPosition;
    } ContactPointState;

    typedef struct {

        std::vector<ContactPointState> leftContactPointsState;
        std::vector<ContactPointState> rightContactPointsState;
        iDynTree::Vector6 momentumInCoM;
        iDynTree::Vector3 comPosition;
        iDynTree::Transform worldToBaseTransform;
        iDynTree::VectorDynSize jointsConfiguration;
        double time;
    } State;

    typedef struct {
        iDynTree::Vector3 pointForceControl;
        iDynTree::Vector3 pointVelocityControl;
    }  ContactPointControl;

    typedef struct {
        std::vector<ContactPointControl> leftContactPointsControl;
        std::vector<ContactPointControl> rightContactPointsControl;
        iDynTree::Vector6 baseVelocityInBaseFrame;
        iDynTree::VectorDynSize jointsVelocity;
        double time;
    } Control;


}

class DynamicalPlanner::Solver{

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    Solver();

    ~Solver();

    bool specifySettings(const Settings& settings);

    bool setInitialCondition(const State& initialState);

    bool setOptimizer(std::shared_ptr<iDynTree::optimization::Optimizer> optimizer);




};

#endif // DPLANNER_SOLVER_H
