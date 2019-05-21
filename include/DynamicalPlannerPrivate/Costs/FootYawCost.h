/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_FOOTYAWCOST_H
#define DPLANNER_FOOTYAWCOST_H

#include <iDynTree/Cost.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <iDynTree/TimeVaryingObject.h>
#include <iDynTree/Core/Position.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class FootYawCost;
    }
}

class DynamicalPlanner::Private::FootYawCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    FootYawCost(const VariablesLabeller& stateVariables, const std::string& footName,
                const std::vector<iDynTree::Position>& pointsPosition);

    ~FootYawCost() override;

    void setDesiredYawTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> desiredYaw);

    virtual bool costEvaluation(double time,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double time,
                                                    const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&,
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

    virtual bool costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure& stateControlSparsity) override;

    virtual bool costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;

};

#endif // DPLANNER_FOOTYAWCOST_H
