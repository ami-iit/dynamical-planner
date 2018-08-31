/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_FRAMEORIENTATIONCOST_H
#define DPLANNER_FRAMEORIENTATIONCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <iDynTree/Core/Rotation.h>
#include <iDynTree/TimeVaryingObject.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class FrameOrientationCost;
    }
}

class DynamicalPlanner::Private::FrameOrientationCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    FrameOrientationCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables, std::shared_ptr<SharedKinDynComputation> sharedKinDyn, const iDynTree::FrameIndex& desiredFrame);

    void setDesiredRotation(const iDynTree::Rotation& desiredRotation);

    bool setDesiredRotationTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingRotation> desiredRotationTrajectory);

    virtual bool costEvaluation(double time,
                                const iDynTree::VectorDynSize& state,
                                const iDynTree::VectorDynSize& control,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double time,
                                                    const iDynTree::VectorDynSize& state,
                                                    const iDynTree::VectorDynSize& control,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&,
                                                      iDynTree::VectorDynSize& partialDerivative) override;
};



#endif // DPLANNER_FRAMEORIENTATIONCOST_H
