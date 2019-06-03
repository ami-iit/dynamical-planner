/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_JOINTSVELOCITYFORPOSTURALCOST_H
#define DPLANNER_JOINTSVELOCITYFORPOSTURALCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <iDynTree/TimeVaryingObject.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class JointsVelocityForPosturalCost;
    }
}

class DynamicalPlanner::Private::JointsVelocityForPosturalCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    JointsVelocityForPosturalCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                                  const iDynTree::VectorDynSize& jointsWeights, const iDynTree::VectorDynSize& posturalGains,
                                  std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredJoints);

    ~JointsVelocityForPosturalCost() override;

    virtual bool costEvaluation(double time,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&control,
                                double& costValue) override;

    virtual bool costFirstPartialDerivativeWRTState(double time,
                                                    const iDynTree::VectorDynSize& state,
                                                    const iDynTree::VectorDynSize& control,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double time, const iDynTree::VectorDynSize& state,
                                                      const iDynTree::VectorDynSize& control,
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

#endif // DPLANNER_JOINTSVELOCITYFORPOSTURALCOST_H
