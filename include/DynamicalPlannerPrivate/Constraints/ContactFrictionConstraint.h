/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CONTACTFRICTIONCONSTRAINT_H
#define DPLANNER_CONTACTFRICTIONCONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/SparsityStructure.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <iDynTree/Core/Position.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        class ContactFrictionConstraint;
    }
}

class DynamicalPlanner::Private::ContactFrictionConstraint : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    ContactFrictionConstraint(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                              const std::string &footName, size_t contactIndex);

    ~ContactFrictionConstraint() override;

    bool setFrictionCoefficient(double frictionCoefficient);

    virtual bool evaluateConstraint(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                    iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                            iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&,
                                              iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;

    virtual bool constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;
};

#endif // DPLANNER_CONTACTFRICTIONCONSTRAINT_H
