/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CONTACTPOSITIONCONSISTENCYCONSTRAINT_H
#define DPLANNER_CONTACTPOSITIONCONSISTENCYCONSTRAINT_H

#include <iDynTree/Constraint.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <private/VariablesLabeller.h>
#include <private/SharedKinDynComputations.h>
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
                                         std::shared_ptr<SharedKinDynComputation> sharedKinDyn,
                                         const iDynTree::FrameIndex &footFrame, const std::string &footName,
                                         const iDynTree::Position &positionInFoot, size_t contactIndex);

    virtual bool evaluateConstraint(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&, iDynTree::MatrixDynSize& jacobian) override;

    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;
};

#endif // DPLANNER_CONTACTPOSITIONCONSISTENCYCONSTRAINT_H
