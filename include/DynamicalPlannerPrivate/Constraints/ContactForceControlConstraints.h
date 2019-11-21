/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CONTACTFORCECONTROLCONSTRAINTS_H
#define DPLANNER_CONTACTFORCECONTROLCONSTRAINTS_H

#include <iDynTree/Constraint.h>
#include <iDynTree/SparsityStructure.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/Utilities/HyperbolicSecant.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        class ContactForceControlConstraints;
    }
}

class DynamicalPlanner::Private::ContactForceControlConstraints : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    ContactForceControlConstraints(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                                   const std::string &footName, size_t contactIndex, const HyperbolicSecant& forceActivation,
                                   double maximumNormalDerivative, double dissipationRatio);

    ~ContactForceControlConstraints() override;

    virtual bool evaluateConstraint(double,
                                    const iDynTree::VectorDynSize& state,
                                    const iDynTree::VectorDynSize& control,
                                    iDynTree::VectorDynSize& constraint) override;

    virtual bool constraintJacobianWRTState(double,
                                            const iDynTree::VectorDynSize& state,
                                            const iDynTree::VectorDynSize& control,
                                            iDynTree::MatrixDynSize& jacobian) override;

    virtual bool constraintJacobianWRTControl(double,
                                              const iDynTree::VectorDynSize& state,
                                              const iDynTree::VectorDynSize& control,
                                              iDynTree::MatrixDynSize& jacobian) override;


    virtual size_t expectedStateSpaceSize() const override;

    virtual size_t expectedControlSpaceSize() const override;

    virtual bool constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;

    virtual bool constraintSecondPartialDerivativeWRTState(double,
                                                           const iDynTree::VectorDynSize& state,
                                                           const iDynTree::VectorDynSize& control,
                                                           const iDynTree::VectorDynSize& lambda,
                                                           iDynTree::MatrixDynSize& hessian) override;

    virtual bool constraintSecondPartialDerivativeWRTControl(double time,
                                                             const iDynTree::VectorDynSize& state,
                                                             const iDynTree::VectorDynSize& control,
                                                             const iDynTree::VectorDynSize& lambda,
                                                             iDynTree::MatrixDynSize& hessian) override;

    virtual bool constraintSecondPartialDerivativeWRTStateControl(double time,
                                                                  const iDynTree::VectorDynSize& state,
                                                                  const iDynTree::VectorDynSize& control,
                                                                  const iDynTree::VectorDynSize& lambda,
                                                                  iDynTree::MatrixDynSize& hessian) override;

    virtual bool constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure& stateSparsity) override;

    virtual bool constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure& stateControlSparsity) override;

    virtual bool constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure& controlSparsity) override;

};


#endif // DPLANNER_CONTACTFORCECONTROLCONSTRAINTS_H
