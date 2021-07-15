/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_STATE_H
#define DPLANNER_STATE_H

#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/TimeVaryingObject.h>
#include <iDynTree/Core/GeomVector3.h>
#include <vector>

namespace DynamicalPlanner {
    class State;

    typedef iDynTree::optimalcontrol::TimeVaryingObject<State> TimeVaryingState;
    typedef iDynTree::optimalcontrol::TimeInvariantObject<DynamicalPlanner::State> TimeInvariantState;

    typedef struct {
        iDynTree::Force pointForce;
        iDynTree::Vector3 pointPosition;
    } ContactPointState;
}

class DynamicalPlanner::State {
public:

    std::vector<ContactPointState> leftContactPointsState;
    std::vector<ContactPointState> rightContactPointsState;
    iDynTree::Vector6 momentumInCoM;
    iDynTree::Vector3 comPosition;
    iDynTree::Transform worldToBaseTransform;
    iDynTree::VectorDynSize jointsConfiguration;
    double time;

    State();

    State(size_t numberOfDofs, size_t numberOfPoints);

    ~State();

    void resize(size_t numberOfDofs, size_t numberOfPoints);

    void zero();

    bool sameSize(const State& other) const;

    bool checkSize(size_t numberOfDofs, size_t numberOfPoints) const;

    iDynTree::Vector3 computeFeetCentroid() const;

    bool leftIsForward(const iDynTree::Vector3& forward);

    double minimumPointForceOnForwardFoot(const iDynTree::Vector3& forwardDirection);

};

extern template class iDynTree::optimalcontrol::TimeInvariantObject<DynamicalPlanner::State>;


#endif // DPLANNER_STATE_H
