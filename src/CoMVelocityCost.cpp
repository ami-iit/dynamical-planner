/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Costs/CoMVelocityCost.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iostream>
#include <cassert>

using namespace DynamicalPlanner::Private;

class LinearMomentumReference : public iDynTree::optimalcontrol::TimeVaryingVector {
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> m_desiredCoMVelocityTrajectory;
    double m_totalMass;
    iDynTree::VectorDynSize m_desiredMomentum;
public:

    LinearMomentumReference(double totalMass) {
        m_totalMass = totalMass;
        m_desiredMomentum.resize(3);
        m_desiredMomentum.zero();
        m_desiredCoMVelocityTrajectory = std::make_shared<iDynTree::optimalcontrol::TimeInvariantVector>(m_desiredMomentum);
    }

    ~LinearMomentumReference() override;

    void setLinearVelocityReference(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMVelocity) {
        m_desiredCoMVelocityTrajectory = desiredCoMVelocity;
    }


    const iDynTree::VectorDynSize& get(double time, bool &isValid) override {

        const iDynTree::VectorDynSize& desiredVelocity = m_desiredCoMVelocityTrajectory->get(time, isValid);

        if (!isValid) {
            m_desiredMomentum.zero();
            return m_desiredMomentum;
        }

        if (desiredVelocity.size() != 3) {
            std::cerr << "[ERROR][LinearMomentumReference::get] The desired CoM velocity at time " << time << " has dimension different from 3.";
            std::cerr << std::endl;
            isValid = false;
            m_desiredMomentum.zero();
            return m_desiredMomentum;
        }

        iDynTree::toEigen(m_desiredMomentum) = m_totalMass * iDynTree::toEigen(desiredVelocity);

        return m_desiredMomentum;
    }
};
LinearMomentumReference::~LinearMomentumReference(){}

class CoMVelocityCost::Implementation{
public:
    std::shared_ptr<LinearMomentumReference> linearMomentumReference;

};

CoMVelocityCost::CoMVelocityCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, unsigned int linearMomentumOffset, std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn)
    : iDynTree::optimalcontrol::L2NormCost ("CoMVelocity", {linearMomentumOffset, 3},
                                            static_cast<unsigned int>(stateVariables.size()),
                                            iDynTree::IndexRange::InvalidRange(),
                                            static_cast<unsigned int>(controlVariables.size()))
    , m_pimpl(std::make_unique<Implementation>())
{
    double totalMass = 0.0;

    const iDynTree::Model & model = timelySharedKinDyn->model();

    for(size_t l=0; l < model.getNrOfLinks(); l++)
    {
        totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
    }

    m_pimpl->linearMomentumReference = std::make_shared<LinearMomentumReference>(totalMass);

    this->setStateDesiredTrajectory(m_pimpl->linearMomentumReference);
}

CoMVelocityCost::~CoMVelocityCost()
{

}

void CoMVelocityCost::setLinearVelocityReference(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMVelocity)
{
    assert(desiredCoMVelocity);
    m_pimpl->linearMomentumReference->setLinearVelocityReference(desiredCoMVelocity);
}
