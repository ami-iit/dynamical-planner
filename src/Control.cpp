/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/Control.h>

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

        template class TimeInvariantObject<Control>;
    }
}


Control::Control()
{ }

Control::Control(size_t numberOfDofs, size_t numberOfPoints)
{
    resize(numberOfDofs, numberOfPoints);
    zero();
}

Control::~Control()
{ }

void Control::resize(size_t numberOfDofs, size_t numberOfPoints)
{
    leftContactPointsControl.resize(numberOfPoints);
    rightContactPointsControl.resize(numberOfPoints);
    jointsVelocity.resize(static_cast<unsigned int>(numberOfDofs));
}

void Control::zero()
{
    for (auto leftPoints : leftContactPointsControl) {
        leftPoints.pointForceControl.zero();
        leftPoints.pointVelocityControl.zero();
    }

    for (auto rightPoints : rightContactPointsControl) {
        rightPoints.pointForceControl.zero();
        rightPoints.pointVelocityControl.zero();
    }

    baseLinearVelocity.zero();
    baseQuaternionDerivative.zero();

    jointsVelocity.zero();

    time = 0.0;
}

bool Control::sameSize(const DynamicalPlanner::Control &other) const
{
    return ((leftContactPointsControl.size() == other.leftContactPointsControl.size())
            && (rightContactPointsControl.size() == other.rightContactPointsControl.size())
            && (jointsVelocity.size() == other.jointsVelocity.size()));
}

bool Control::checkSize(size_t numberOfDofs, size_t numberOfPoints) const
{
    return ((leftContactPointsControl.size() == numberOfPoints)
            && (rightContactPointsControl.size() == numberOfPoints)
            && (jointsVelocity.size() == numberOfDofs));
}
