/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SWINGCOST_H
#define DPLANNER_SWINGCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class SwingCost;
    }
}

class DynamicalPlanner::Private::SwingCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    SwingCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
              const std::string &footName, size_t contactIndex, double desiredSwingHeight);

    ~SwingCost() override;

    virtual bool costEvaluation(double,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&control,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double,
                                                    const iDynTree::VectorDynSize& state,
                                                    const iDynTree::VectorDynSize& control,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize&state, const iDynTree::VectorDynSize&control,
                                                      iDynTree::VectorDynSize& partialDerivative) override;
};

#endif // DPLANNER_SWINGCOST_H
