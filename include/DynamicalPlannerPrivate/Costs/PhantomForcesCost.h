/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_PHANTOMFORCESCOST_H
#define DPLANNER_PHANTOMFORCESCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/Utilities/HyperbolicSecant.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class PhantomForcesCost;
    }
}

class DynamicalPlanner::Private::PhantomForcesCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    PhantomForcesCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                      const std::string &footName, size_t contactIndex, const HyperbolicSecant& forceActivation);

    ~PhantomForcesCost() override;

    virtual bool costEvaluation(double,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&control,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double,
                                                    const iDynTree::VectorDynSize& state,
                                                    const iDynTree::VectorDynSize& control,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&,
                                                      iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costSecondPartialDerivativeWRTState(double time,
                                                     const iDynTree::VectorDynSize& state,
                                                     const iDynTree::VectorDynSize& control,
                                                     iDynTree::MatrixDynSize& partialDerivative) override;

    virtual bool costSecondPartialDerivativeWRTControl(double time,
                                                       const iDynTree::VectorDynSize& state,
                                                       const iDynTree::VectorDynSize& control,
                                                       iDynTree::MatrixDynSize& partialDerivative) override;

    virtual bool costSecondPartialDerivativeWRTStateControl(double time,
                                                            const iDynTree::VectorDynSize& state,
                                                            const iDynTree::VectorDynSize& control,
                                                            iDynTree::MatrixDynSize& partialDerivative) override;
};

#endif // DPLANNER_PHANTOMFORCESCOST_H
