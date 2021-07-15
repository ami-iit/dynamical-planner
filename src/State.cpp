/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/State.h>
#include <iDynTree/Core/EigenHelpers.h>

using namespace DynamicalPlanner;


namespace iDynTree {
    namespace optimalcontrol {

        template<typename Object>
        TimeInvariantObject<Object>::TimeInvariantObject()
        { }



        template<typename Object>
        TimeInvariantObject<Object>::TimeInvariantObject(const Object &timeInvariantObject)
            : m_invariantObject(timeInvariantObject)
        { }

        template<typename Object>
        TimeInvariantObject<Object>::~TimeInvariantObject()
        { }

        template<typename Object>
        Object &TimeInvariantObject<Object>::get()
        {
            return m_invariantObject;
        }

        template<typename Object>
        const Object &TimeInvariantObject<Object>::get(double /*time*/, bool &isValid)
        {
            isValid = true;
            return m_invariantObject;
        }

        template class TimeInvariantObject<State>;
    }
}


State::State()
{ }

State::State(size_t numberOfDofs, size_t numberOfPoints)
{
    resize(numberOfDofs, numberOfPoints);
    zero();
}

State::~State()
{ }

void State::resize(size_t numberOfDofs, size_t numberOfPoints)
{
    leftContactPointsState.resize(numberOfPoints);
    rightContactPointsState.resize(numberOfPoints);
    jointsConfiguration.resize(static_cast<unsigned int>(numberOfDofs));
}

void State::zero()
{
    for (auto& leftPoint : leftContactPointsState) {
        leftPoint.pointForce.zero();
        leftPoint.pointPosition.zero();
    }

    for (auto& rightPoint : rightContactPointsState) {
        rightPoint.pointForce.zero();
        rightPoint.pointPosition.zero();
    }

    momentumInCoM.zero();
    comPosition.zero();
    worldToBaseTransform = iDynTree::Transform::Identity();
    jointsConfiguration.zero();
    time = 0;
}

bool State::sameSize(const DynamicalPlanner::State &other) const
{
    return ((leftContactPointsState.size() == other.leftContactPointsState.size())
            && (rightContactPointsState.size() == other.rightContactPointsState.size())
            && (jointsConfiguration.size() == other.jointsConfiguration.size()));
}

bool State::checkSize(size_t numberOfDofs, size_t numberOfPoints) const
{
    return ((leftContactPointsState.size() == numberOfPoints)
            && (rightContactPointsState.size() == numberOfPoints)
            && (jointsConfiguration.size() == numberOfDofs));
}

iDynTree::Vector3 State::computeFeetCentroid() const{
    iDynTree::Vector3 meanPosition;
    meanPosition.zero();

    if (leftContactPointsState.size() + rightContactPointsState.size() == 0)
    {
        return meanPosition;
    }

    for (size_t i = 0; i < leftContactPointsState.size(); ++i) {
        iDynTree::toEigen(meanPosition) += iDynTree::toEigen(leftContactPointsState[i].pointPosition);
    }

    for (size_t i = 0; i < rightContactPointsState.size(); ++i) {
        iDynTree::toEigen(meanPosition) += iDynTree::toEigen(rightContactPointsState[i].pointPosition);
    }

    iDynTree::toEigen(meanPosition) /= leftContactPointsState.size() + rightContactPointsState.size();

    return meanPosition;
}
