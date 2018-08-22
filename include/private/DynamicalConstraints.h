/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_DYNAMICALCONSTRAINTS_H
#define DPLANNER_DYNAMICALCONSTRAINTS_H

#include <iDynTree/DynamicalSystem.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/Position.h>
#include <private/VariablesLabeller.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class DynamicalConstraints;
    }
}

class DynamicalPlanner::Private::DynamicalConstraints : public iDynTree::optimalcontrol::DynamicalSystem {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    DynamicalConstraints(VariablesLabeller& stateVariables, VariablesLabeller& controlVariables, double totalMass);

    ~DynamicalConstraints() override;

    virtual bool dynamics(const iDynTree::VectorDynSize& state, double,
                          iDynTree::VectorDynSize& stateDynamics) override;

    virtual bool dynamicsStateFirstDerivative(const iDynTree::VectorDynSize& state, double,
                                              iDynTree::MatrixDynSize& dynamicsDerivative) override;

    virtual bool dynamicsControlFirstDerivative(const iDynTree::VectorDynSize& state,
                                                double time,
                                                iDynTree::MatrixDynSize& dynamicsDerivative) override;

};

#endif // DPLANNER_DYNAMICALCONSTRAINTS_H
