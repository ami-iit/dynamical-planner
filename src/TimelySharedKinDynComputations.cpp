/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/TimelySharedKinDynComputations.h>
#include <iDynTree/Core/Utils.h>
#include <cmath>
#include <cassert>

using namespace DynamicalPlanner::Private;



size_t TimelySharedKinDynComputations::getClosestIndex(double time)
{
    // Corner cases
    if (time <= m_pointerContainer[0].time) {
        m_previousIndex = 0;
        return 0;
    }
    if (time >= m_pointerContainer.back().time) {
        m_previousIndex = m_pointerContainer.size() - 1;
        return m_pointerContainer.size() - 1;
    }

    if (iDynTree::checkDoublesAreEqual(time, m_pointerContainer[m_previousIndex].time)) {
        return m_previousIndex;
    }

    if (iDynTree::checkDoublesAreEqual(time, m_pointerContainer[m_previousIndex + 1].time)) {
        m_previousIndex++;
        return m_previousIndex;
    }

    // Doing binary search
    size_t i = 0, j = m_pointerContainer.size(), mid = 0;
    while (i < j) {
        mid = (i + j) / 2;

        if (iDynTree::checkDoublesAreEqual(time, m_pointerContainer[mid].time)) {
            m_previousIndex = mid;
            return mid;
        }

        /* If target is less than array element,
             then search in left */
        if (time < m_pointerContainer[mid].time) {

            // If target is greater than previous
            // to mid, return closest of two
            if (time > m_pointerContainer[mid - 1].time) {
                m_previousIndex =  ((time - m_pointerContainer[mid - 1].time) >= (m_pointerContainer[mid].time - time)) ? mid : mid - 1;
                return m_previousIndex;
            }

            /* Repeat for left half */
            j = mid;
        }

        // If target is greater than mid
        else {
            if (time < m_pointerContainer[mid + 1].time) {
                m_previousIndex = ((time - m_pointerContainer[mid].time) >= (m_pointerContainer[mid + 1].time - time)) ? mid + 1 : mid;
                return m_previousIndex;
            }
            // update i
            i = mid + 1;
        }
    }

    // Only single element left after search
    m_previousIndex = mid;
    return mid;
}

TimelySharedKinDynComputations::TimelySharedKinDynComputations()
    : m_timingsSet(false)
    , m_previousIndex(0)
{
    m_sharedTemplate = std::make_shared<SharedKinDynComputations>();
}

TimelySharedKinDynComputations::~TimelySharedKinDynComputations()
{ }

bool TimelySharedKinDynComputations::loadRobotModel(const iDynTree::Model &model)
{
    return m_sharedTemplate->loadRobotModel(model);
}

const iDynTree::Model &TimelySharedKinDynComputations::model() const
{
    return m_sharedTemplate->model();
}

bool TimelySharedKinDynComputations::isValid() const
{
    return m_sharedTemplate->isValid();
}

void TimelySharedKinDynComputations::setGravity(const iDynTree::Vector3 &gravity)
{
    return m_sharedTemplate->setGravity(gravity);
}

const iDynTree::Vector3 &TimelySharedKinDynComputations::gravity() const
{
    return m_sharedTemplate->gravity();
}

bool TimelySharedKinDynComputations::setToleranceForUpdate(double tol)
{
    return m_sharedTemplate->setToleranceForUpdate(tol);
}

double TimelySharedKinDynComputations::getUpdateTolerance() const
{
    return m_sharedTemplate->getUpdateTolerance();
}

bool TimelySharedKinDynComputations::setFloatingBase(const std::string &floatingBaseName)
{
    return m_sharedTemplate->setFloatingBase(floatingBaseName);
}

bool TimelySharedKinDynComputations::setTimings(const std::vector<double> &timings)
{
    if (!isValid()) {
        return false;
    }

    for (size_t i = 0; i < timings.size(); ++i) {
        if (i < m_pointerContainer.size()) {
            m_pointerContainer[i].time = timings[i];
            //the pointer is not updated
        } else {
            TimedSharedKinDyn newPointer;
            newPointer.time = timings[i];
            newPointer.pointer = std::make_shared<SharedKinDynComputations>(*m_sharedTemplate);
            m_pointerContainer.push_back(newPointer);
        }
    }

    m_timingsSet = true;
    return true;
}

SharedKinDynComputationsPointer TimelySharedKinDynComputations::get(double time)
{
    assert(isValid());

    if (!m_timingsSet) {
        return m_sharedTemplate;
    }

    return m_pointerContainer[getClosestIndex(time)].pointer;
}
