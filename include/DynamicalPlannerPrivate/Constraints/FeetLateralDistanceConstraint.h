/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_FEETLATERALDISTANCECONSTRAINT_H
#define DPLANNER_FEETLATERALDISTANCECONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/SparsityStructure.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/TimelySharedKinDynComputations.h>
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
                                  std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, unsigned int lateralIndex,
                                  iDynTree::FrameIndex referenceFootFrame, iDynTree::FrameIndex otherFootFrame);

    ~FeetLateralDistanceConstraint() override;

    bool setMinimumDistance(double minDistance);

    virtual bool evaluateConstraint(double time, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double time, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;

    virtual bool constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;
};



#endif // DPLANNER_FEETLATERALDISTANCECONSTRAINT_H
