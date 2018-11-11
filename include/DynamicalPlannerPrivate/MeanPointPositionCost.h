/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_MEANPOINTPOSITIONCOST_H
#define DPLANNER_MEANPOINTPOSITIONCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <iDynTree/TimeVaryingObject.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class MeanPointPositionCost;
    }
}

class DynamicalPlanner::Private::MeanPointPositionCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    MeanPointPositionCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables);

    ~MeanPointPositionCost() override;

    bool setDesiredPositionTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> desiredPosition);

    virtual bool costEvaluation(double time,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double time,
                                                    const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&,
                                                      iDynTree::VectorDynSize& partialDerivative) override;
};

#endif // DPLANNER_MEANPOINTPOSITIONCOST_H
