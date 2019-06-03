/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/JointsVelocityForPosturalCost.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class JointsVelocityForPosturalCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;
    iDynTree::IndexRange jointsVelocityRange, jointsRange;

    iDynTree::VectorDynSize posturalGains, jointsBuffer, jointsError, posturalGainsSquared,
        jointsWeights, jointsWeightsSquared, weightsTimesGains;
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredJoints;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};

JointsVelocityForPosturalCost::JointsVelocityForPosturalCost(const VariablesLabeller &stateVariables,
                                                             const VariablesLabeller &controlVariables,
                                                             const iDynTree::VectorDynSize &jointsWeights,
                                                             const iDynTree::VectorDynSize &posturalGains,
                                                             std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredJoints)
    : iDynTree::optimalcontrol::Cost ("JointsVelocityForPostural")
      , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->jointsRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsRange.isValid());

    m_pimpl->jointsVelocityRange = m_pimpl->controlVariables.getIndexRange("JointsVelocity");
    assert(m_pimpl->jointsVelocityRange.isValid());

    m_pimpl->posturalGains = posturalGains;
    assert(m_pimpl->posturalGains.size() == m_pimpl->jointsRange.size);

    m_pimpl->jointsWeights = jointsWeights;
    assert(m_pimpl->jointsWeights.size() == m_pimpl->jointsRange.size);

    for (unsigned int i = 0; i < m_pimpl->posturalGains.size(); ++i) {
        m_pimpl->posturalGains(i) = m_pimpl->jointsWeights(i) * m_pimpl->posturalGains(i);
    }

    m_pimpl->desiredJoints = desiredJoints;
    assert(m_pimpl->desiredJoints);

    m_pimpl->jointsBuffer.resize(static_cast<unsigned int>(m_pimpl->jointsRange.size));
    m_pimpl->jointsError.resize(static_cast<unsigned int>(m_pimpl->jointsRange.size));
    m_pimpl->posturalGainsSquared.resize(static_cast<unsigned int>(m_pimpl->jointsRange.size));
    m_pimpl->jointsWeightsSquared.resize(static_cast<unsigned int>(m_pimpl->jointsRange.size));
    m_pimpl->weightsTimesGains.resize(static_cast<unsigned int>(m_pimpl->jointsRange.size));

    for (unsigned int i = 0; i < m_pimpl->posturalGains.size(); ++i) {
        m_pimpl->posturalGainsSquared(i) = m_pimpl->posturalGains(i) * m_pimpl->posturalGains(i);
        m_pimpl->jointsWeightsSquared(i) = m_pimpl->jointsWeights(i) * m_pimpl->jointsWeights(i);
        m_pimpl->weightsTimesGains(i) = m_pimpl->posturalGains(i) * m_pimpl->jointsWeights(i);
    }

    m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->jointsRange.offset),
                                                   static_cast<size_t>(m_pimpl->jointsRange.offset),
                                                   static_cast<size_t>(m_pimpl->jointsRange.size));

    m_pimpl->mixedHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->jointsRange.offset),
                                                   static_cast<size_t>(m_pimpl->jointsVelocityRange.offset),
                                                   static_cast<size_t>(m_pimpl->jointsRange.size));

    m_pimpl->controlHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->jointsVelocityRange.offset),
                                                     static_cast<size_t>(m_pimpl->jointsVelocityRange.offset),
                                                     static_cast<size_t>(m_pimpl->jointsVelocityRange.size));
}

JointsVelocityForPosturalCost::~JointsVelocityForPosturalCost()
{

}

bool JointsVelocityForPosturalCost::costEvaluation(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    bool isValid = false;
    const iDynTree::VectorDynSize& desiredJoints = m_pimpl->desiredJoints->get(time, isValid);
    assert(isValid);

    iDynTree::toEigen(m_pimpl->jointsError) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->jointsRange)) -
        iDynTree::toEigen(desiredJoints);

    iDynTree::toEigen(m_pimpl->jointsBuffer) = iDynTree::toEigen(m_pimpl->jointsWeights).asDiagonal() *
            iDynTree::toEigen(m_pimpl->controlVariables(m_pimpl->jointsVelocityRange)) +
        iDynTree::toEigen(m_pimpl->posturalGains).asDiagonal() * iDynTree::toEigen(m_pimpl->jointsError);

    costValue = 0.5 * iDynTree::toEigen(m_pimpl->jointsBuffer).squaredNorm();

    return true;
}

bool JointsVelocityForPosturalCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    bool isValid = false;
    const iDynTree::VectorDynSize& desiredJoints = m_pimpl->desiredJoints->get(time, isValid);
    assert(isValid);

    iDynTree::toEigen(m_pimpl->jointsError) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->jointsRange)) -
        iDynTree::toEigen(desiredJoints);

    iDynTree::toEigen(m_pimpl->jointsBuffer) = iDynTree::toEigen(m_pimpl->jointsWeights).asDiagonal() *
            iDynTree::toEigen(m_pimpl->controlVariables(m_pimpl->jointsVelocityRange)) +
        iDynTree::toEigen(m_pimpl->posturalGains).asDiagonal() * iDynTree::toEigen(m_pimpl->jointsError);

    iDynTree::toEigen(partialDerivative).segment(m_pimpl->jointsRange.offset, m_pimpl->jointsRange.size) =
        iDynTree::toEigen(m_pimpl->posturalGains).asDiagonal() * iDynTree::toEigen(m_pimpl->jointsBuffer);

    return true;
}

bool JointsVelocityForPosturalCost::costFirstPartialDerivativeWRTControl(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    bool isValid = false;
    const iDynTree::VectorDynSize& desiredJoints = m_pimpl->desiredJoints->get(time, isValid);
    assert(isValid);

    iDynTree::toEigen(m_pimpl->jointsError) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->jointsRange)) -
        iDynTree::toEigen(desiredJoints);

    iDynTree::toEigen(m_pimpl->jointsBuffer) = iDynTree::toEigen(m_pimpl->jointsWeights).asDiagonal() *
            iDynTree::toEigen(m_pimpl->controlVariables(m_pimpl->jointsVelocityRange)) +
        iDynTree::toEigen(m_pimpl->posturalGains).asDiagonal() * iDynTree::toEigen(m_pimpl->jointsError);

    iDynTree::toEigen(partialDerivative).segment(m_pimpl->jointsVelocityRange.offset, m_pimpl->jointsVelocityRange.size) =
        iDynTree::toEigen(m_pimpl->jointsWeights).asDiagonal() * iDynTree::toEigen(m_pimpl->jointsBuffer);

    return true;
}

bool JointsVelocityForPosturalCost::costSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    iDynTree::toEigen(partialDerivative).block(m_pimpl->jointsRange.offset, m_pimpl->jointsRange.offset,
                                               m_pimpl->jointsRange.size, m_pimpl->jointsRange.size) =
        iDynTree::toEigen(m_pimpl->posturalGainsSquared).asDiagonal();

    return true;
}

bool JointsVelocityForPosturalCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    iDynTree::toEigen(partialDerivative).block(m_pimpl->jointsVelocityRange.offset, m_pimpl->jointsVelocityRange.offset,
                                               m_pimpl->jointsVelocityRange.size, m_pimpl->jointsVelocityRange.size) =
        iDynTree::toEigen(m_pimpl->jointsWeightsSquared).asDiagonal();
    return true;
}

bool JointsVelocityForPosturalCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    iDynTree::toEigen(partialDerivative).block(m_pimpl->jointsRange.offset, m_pimpl->jointsVelocityRange.offset,
                                               m_pimpl->jointsRange.size, m_pimpl->jointsVelocityRange.size) =
        iDynTree::toEigen(m_pimpl->weightsTimesGains).asDiagonal();

    return true;
}

bool JointsVelocityForPosturalCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool JointsVelocityForPosturalCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool JointsVelocityForPosturalCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
