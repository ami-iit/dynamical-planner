/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_PLANARVELOCITYCONTROLCONSTRAINTS_H
#define DPLANNER_PLANARVELOCITYCONTROLCONSTRAINTS_H

#include <iDynTree/Constraint.h>
#include <iDynTree/SparsityStructure.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/HyperbolicTangent.h>
#include <DynamicalPlannerPrivate/HyperbolicSecant.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        class PlanarVelocityControlConstraints;
    }
}

class DynamicalPlanner::Private::PlanarVelocityControlConstraints : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    PlanarVelocityControlConstraints(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                                     const std::string &footName, size_t contactIndex, const HyperbolicTangent& planarVelocityActivation,
                                     double xMaximumDerivative, double yMaximumDerivative);

    ~PlanarVelocityControlConstraints() override;

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

};

#endif // DPLANNER_PLANARVELOCITYCONTROLCONSTRAINTS_H
