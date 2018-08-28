/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CONTACTVELOCITYCONTROLCONSTRAINTS_H
#define DPLANNER_CONTACTVELOCITYCONTROLCONSTRAINTS_H

#include <iDynTree/Constraint.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/MatrixDynSize.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/HyperbolicSecant.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {
        class ContactVelocityControlConstraints;
    }
}

class DynamicalPlanner::Private::ContactVelocityControlConstraints : public iDynTree::optimalcontrol::Constraint {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    ContactVelocityControlConstraints(const VariablesLabeller& stateVariables, const VariablesLabeller& controlVariables,
                                      const std::string &footName, size_t contactIndex, const HyperbolicSecant& activationXY,
                                      const HyperbolicSecant& activationZ,
                                      const iDynTree::Vector3& maximumDerivatives, const iDynTree::Vector3& dissipationRatios);

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

};


#endif // DPLANNER_CONTACTVELOCITYCONTROLCONSTRAINTS_H