/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CONTACTPOSITIONCONSISTENCYCONSTRAINT_H
#define DPLANNER_CONTACTPOSITIONCONSISTENCYCONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/SparsityStructure.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/TimelySharedKinDynComputations.h>
#include <iDynTree/Core/Position.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        class ContactPositionConsistencyConstraint;
    }
}

class DynamicalPlanner::Private::ContactPositionConsistencyConstraint : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    ContactPositionConsistencyConstraint(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                                         std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                         iDynTree::FrameIndex footFrame, const std::string &footName,
                                         const iDynTree::Position &positionInFoot, size_t contactIndex);

    ~ContactPositionConsistencyConstraint() override;

    void setEqualityTolerance(double tolerance);

    virtual bool evaluateConstraint(double time, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double time, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;

    virtual bool constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;
};

#endif // DPLANNER_CONTACTPOSITIONCONSISTENCYCONSTRAINT_H
