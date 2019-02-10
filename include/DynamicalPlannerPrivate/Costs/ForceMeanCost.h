/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_FORCEMEANCOST_H
#define DPLANNER_FORCEMEANCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class ForceMeanCost;
    }
}

class DynamicalPlanner::Private::ForceMeanCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    ForceMeanCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                  const std::string &footName, size_t contactIndex);

    ~ForceMeanCost() override;

    virtual bool costEvaluation(double,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double,
                                                    const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&,
                                                      iDynTree::VectorDynSize& partialDerivative) override;
};


#endif // DPLANNER_FORCEMEANCOST_H
