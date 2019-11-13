/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_FEETRELATIVEHEIGHTCONSTRAINT_H
#define DPLANNER_FEETRELATIVEHEIGHTCONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/SparsityStructure.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <memory>

namespace DynamicalPlanner {
namespace Private {
class FeetRelativeHeightConstraint;
}
}

class DynamicalPlanner::Private::FeetRelativeHeightConstraint : public  iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    FeetRelativeHeightConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                 double maxRelativeHeight);

    ~FeetRelativeHeightConstraint() override;

    virtual bool evaluateConstraint(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;

    virtual bool constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;

    virtual bool constraintSecondPartialDerivativeWRTState(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize&) override;

    virtual bool constraintSecondPartialDerivativeWRTControl(double time,
                                                             const iDynTree::VectorDynSize& state,
                                                             const iDynTree::VectorDynSize& control,
                                                             const iDynTree::VectorDynSize& lambda,
                                                             iDynTree::MatrixDynSize& hessian) override;

    virtual bool constraintSecondPartialDerivativeWRTStateControl(double time,
                                                                  const iDynTree::VectorDynSize& state,
                                                                  const iDynTree::VectorDynSize& control,
                                                                  const iDynTree::VectorDynSize& lambda,
                                                                  iDynTree::MatrixDynSize& hessian) override;

    virtual bool constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure& stateControlSparsity) override;

    virtual bool constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;
};


#endif // DPLANNER_FEETRELATIVEHEIGHTCONSTRAINT_H
