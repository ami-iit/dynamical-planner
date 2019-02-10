/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_TIMELYSHAREDKINDYNCOMPUTATIONS_H
#define DPLANNER_TIMELYSHAREDKINDYNCOMPUTATIONS_H

#include <DynamicalPlannerPrivate/Utilities/SharedKinDynComputations.h>
#include <memory>
#include <vector>

namespace DynamicalPlanner {
    namespace Private {
        class TimelySharedKinDynComputations;
    }
}

class DynamicalPlanner::Private::TimelySharedKinDynComputations {
    typedef struct {
        SharedKinDynComputationsPointer pointer;
        double time;
    } TimedSharedKinDyn;

    SharedKinDynComputationsPointer m_sharedTemplate;
    std::vector<TimedSharedKinDyn> m_pointerContainer;
    bool m_timingsSet;
    size_t m_previousIndex;

    size_t getClosestIndex(double time);

public:

    TimelySharedKinDynComputations();

    ~TimelySharedKinDynComputations();

    bool loadRobotModel(const iDynTree::Model& model);

    const iDynTree::Model& model() const;

    bool isValid() const;

    void setGravity(const iDynTree::Vector3& gravity);

    const iDynTree::Vector3 &gravity() const;

    bool setToleranceForUpdate(double tol);

    double getUpdateTolerance() const;

    bool setFloatingBase(const std::string & floatingBaseName);

    bool setTimings(const std::vector<double> &timings);

    SharedKinDynComputationsPointer get(double time);
};


#endif // DPLANNER_TIMELYSHAREDKINDYNCOMPUTATIONS_H
