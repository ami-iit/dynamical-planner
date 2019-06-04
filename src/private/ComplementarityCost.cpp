/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/ComplementarityCost.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::ComplementarityCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;

    iDynTree::IndexRange velocityPointRange, forcePointRange;
    unsigned int vzCol, fzCol;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

};



ComplementarityCost::ComplementarityCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                         const std::string &footName, size_t contactIndex)
    : iDynTree::optimalcontrol::Cost ("Complementarity" + footName + std::to_string(contactIndex))
      , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->velocityPointRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->velocityPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->vzCol = static_cast<unsigned int>(m_pimpl->velocityPointRange.offset + 2);
    m_pimpl->fzCol = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    m_pimpl->stateHessianSparsity.add(m_pimpl->fzCol, m_pimpl->fzCol);
    m_pimpl->mixedHessianSparsity.add(m_pimpl->fzCol, m_pimpl->vzCol);
    m_pimpl->controlHessianSparsity.add(m_pimpl->vzCol, m_pimpl->vzCol);
}

ComplementarityCost::~ComplementarityCost()
{ }

bool ComplementarityCost::costEvaluation(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control,
                                         double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    double vz = m_pimpl->controlVariables(m_pimpl->velocityPointRange)(2);
    double fz = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);

    costValue = 0.5 * vz * vz * fz * fz;

    return true;
}

bool ComplementarityCost::costFirstPartialDerivativeWRTState(double, const iDynTree::VectorDynSize &state,
                                                             const iDynTree::VectorDynSize &control,
                                                             iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    double vz = m_pimpl->controlVariables(m_pimpl->velocityPointRange)(2);
    double fz = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);

    partialDerivative(m_pimpl->fzCol) = vz * vz * fz;

    return true;
}

bool ComplementarityCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &state,
                                                               const iDynTree::VectorDynSize &control,
                                                               iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    double vz = m_pimpl->controlVariables(m_pimpl->velocityPointRange)(2);
    double fz = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);

    partialDerivative(m_pimpl->vzCol) = vz * fz * fz;

    return true;
}

bool ComplementarityCost::costSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                              const iDynTree::VectorDynSize &control,
                                                              iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->controlVariables = control;

    double vz = m_pimpl->controlVariables(m_pimpl->velocityPointRange)(2);

    partialDerivative(m_pimpl->fzCol, m_pimpl->fzCol) = vz*vz;

    return true;
}

bool ComplementarityCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &state,
                                                                const iDynTree::VectorDynSize &/*control*/,
                                                                iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    double fz = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);

    partialDerivative (m_pimpl->vzCol, m_pimpl->vzCol) = fz * fz;

    return true;
}

bool ComplementarityCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &state,
                                                                     const iDynTree::VectorDynSize &control,
                                                                     iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    double vz = m_pimpl->controlVariables(m_pimpl->velocityPointRange)(2);
    double fz = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);

    partialDerivative (m_pimpl->fzCol, m_pimpl->vzCol) = 2 * fz * vz;

    return true;
}

bool ComplementarityCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ComplementarityCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ComplementarityCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
