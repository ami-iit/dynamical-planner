/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_STATICTORQUESCOST_H
#define DPLANNER_STATICTORQUESCOST_H

#include <iDynTree/Cost.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class StaticTorquesCost;
    }
}

class DynamicalPlanner::Private::StaticTorquesCost : public iDynTree::optimalcontrol::Cost {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    StaticTorquesCost(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                      std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, const iDynTree::FrameIndex &leftFootFrame,
                      const iDynTree::FrameIndex &rightFootFrame, const std::vector<iDynTree::Position> &positionsInLeftFoot,
                      const std::vector<iDynTree::Position> &positionsInRightFoot);

    ~StaticTorquesCost() override;

    bool setWeights(const iDynTree::VectorDynSize& torquesWeights);

    void computeStaticTorques(double time, const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize& staticTorques); //DEBUG

    virtual bool costEvaluation(double time,
                                const iDynTree::VectorDynSize& state, const iDynTree::VectorDynSize&control,
                                double& costValue) override;

    void computeStaticTorquesJacobian(double time, const iDynTree::VectorDynSize& state,
                                      const iDynTree::VectorDynSize& control,
                                      iDynTree::MatrixDynSize& staticTorquesJacobian); //DEBUG

    virtual bool costFirstPartialDerivativeWRTState(double time,
                                                    const iDynTree::VectorDynSize& state,
                                                    const iDynTree::VectorDynSize& control,
                                                    iDynTree::VectorDynSize& partialDerivative) override;

    virtual bool costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize&, const iDynTree::VectorDynSize&,
                                                      iDynTree::VectorDynSize& partialDerivative) override;
};

#endif // DPLANNER_STATICTORQUESCOST_H
