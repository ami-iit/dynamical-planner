/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_COMVELOCITYCOST_H
#define DPLANNER_COMVELOCITYCOST_H

#include <iDynTree/L2NormCost.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/TimelySharedKinDynComputations.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class CoMVelocityCost;
    }
}

class DynamicalPlanner::Private::CoMVelocityCost : public iDynTree::optimalcontrol::L2NormCost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    CoMVelocityCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                    unsigned int linearMomentumOffset, std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn);

    ~CoMVelocityCost() override;

    void setLinearVelocityReference(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMVelocity);

};


#endif // DPLANNER_COMVELOCITYCOST_H
