/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_FEETLATERALDISTANCECONSTRAINT_H
#define DPLANNER_FEETLATERALDISTANCECONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <private/VariablesLabeller.h>
#include <private/SharedKinDynComputations.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class FeetLateralDistanceConstraint;
    }
}

class DynamicalPlanner::Private::FeetLateralDistanceConstraint : public  iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    FeetLateralDistanceConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                  std::shared_ptr<SharedKinDynComputation> sharedKinDyn, unsigned int lateralIndex,
                                  iDynTree::FrameIndex referenceFootFrame, iDynTree::FrameIndex otherFootFrame);

    bool setMinimumDistance(double minDistance);

    virtual bool evaluateConstraint(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;
};



#endif // DPLANNER_FEETLATERALDISTANCECONSTRAINT_H
