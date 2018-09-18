/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CONTROL_H
#define DPLANNER_CONTROL_H

#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/TimeVaryingObject.h>


namespace DynamicalPlanner {
    class Control;

    typedef iDynTree::optimalcontrol::TimeVaryingObject<Control> TimeVaryingControl;
    typedef iDynTree::optimalcontrol::TimeInvariantObject<DynamicalPlanner::Control> TimeInvariantControl;


    typedef struct {
        iDynTree::Vector3 pointForceControl;
        iDynTree::Vector3 pointVelocityControl;
    }  ContactPointControl;
}

class DynamicalPlanner::Control {
public:

    std::vector<ContactPointControl> leftContactPointsControl;
    std::vector<ContactPointControl> rightContactPointsControl;
    iDynTree::Vector6 baseVelocityInBaseFrame;
    iDynTree::VectorDynSize jointsVelocity;
    double time;

    Control();

    Control(size_t numberOfDofs, size_t numberOfPoints);

    ~Control();

    void resize(size_t numberOfDofs, size_t numberOfPoints);

    void zero();

    bool sameSize(const Control& other) const;

    bool checkSize(size_t numberOfDofs, size_t numberOfPoints) const;

};

extern template class iDynTree::optimalcontrol::TimeInvariantObject<DynamicalPlanner::Control>;

#endif // DPLANNER_CONTROL_H
