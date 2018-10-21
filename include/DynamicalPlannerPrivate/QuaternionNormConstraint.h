/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_QUATERNIONNORMCONSTRAINT_H
#define DPLANNER_QUATERNIONNORMCONSTRAINT_H

#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <iDynTree/SparsityStructure.h>
#include <iDynTree/Constraint.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class QuaternionNormConstraint;
    }
}

class DynamicalPlanner::Private::QuaternionNormConstraint : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    QuaternionNormConstraint(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables);

    ~QuaternionNormConstraint() override;

    void setEqualityTolerance(double tolerance);

    virtual bool evaluateConstraint(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;

    virtual bool constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;
};

#endif // DPLANNER_QUATERNIONNORMCONSTRAINT_H
