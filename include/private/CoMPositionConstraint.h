/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_COMPOSITIONCONSTRAINT_H
#define DPLANNER_COMPOSITIONCONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/Model/Model.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <private/VariablesLabeller.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        class CoMPositionConstraint;
    }
}

class DynamicalPlanner::Private::CoMPositionConstraint : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    CoMPositionConstraint(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables, const iDynTree::Model &model, const std::string &floatingBase);

    virtual bool evaluateConstraint(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;
};

#endif // DPLANNER_COMPOSITIONCONSTRAINT_H