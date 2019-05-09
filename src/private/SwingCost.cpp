/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/SwingCost.h>

using namespace DynamicalPlanner::Private;

class SwingCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;
    iDynTree::IndexRange pointVelocityRange, pointPositionRange;
    double swingHeight;

    std::string footName;
    size_t contactIndex;

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};

SwingCost::SwingCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, const std::string &footName, size_t contactIndex, double desiredSwingHeight)
    : iDynTree::optimalcontrol::Cost ("Swing" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->pointVelocityRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->pointVelocityRange.isValid());
    m_pimpl->pointPositionRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->pointPositionRange.isValid());

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();

    m_pimpl->swingHeight = desiredSwingHeight;

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);

    m_pimpl->stateHessianSparsity.add(pzIndex, pzIndex);
    m_pimpl->controlHessianSparsity.add(uxIndex, uxIndex);
    m_pimpl->controlHessianSparsity.add(uyIndex, uyIndex);
    m_pimpl->mixedHessianSparsity.add(pzIndex, uxIndex);
    m_pimpl->mixedHessianSparsity.add(pzIndex, uyIndex);

}

SwingCost::~SwingCost()
{ }

bool SwingCost::costEvaluation(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    double costX = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);
    double costY = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    costValue = 0.5 * (costX * costX + costY * costY);
    return true;
}

bool SwingCost::costFirstPartialDerivativeWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    double costX = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);
    double costY = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    m_pimpl->stateGradientBuffer(pzIndex) = costX * m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) + costY * m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1);

    partialDerivative = m_pimpl->stateGradientBuffer;

    return true;
}

bool SwingCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);

    double costX = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);
    double costY = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    m_pimpl->controlGradientBuffer(uxIndex) = costX * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);
    m_pimpl->controlGradientBuffer(uyIndex) = costY * (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    partialDerivative = m_pimpl->controlGradientBuffer;

    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    double ux = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0);
    double uy = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1);

    partialDerivative(pzIndex, pzIndex) = ux * ux + uy * uy;

    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);
    double heightDifference = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    partialDerivative(uxIndex, uxIndex) = heightDifference * heightDifference;
    partialDerivative(uyIndex, uyIndex) = partialDerivative(uxIndex, uxIndex);
    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);
    double ux = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0);
    double uy = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1);
    double heightDifference = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    partialDerivative(pzIndex, uxIndex) = 2.0 * ux * heightDifference;
    partialDerivative(pzIndex, uyIndex) = 2.0 * uy * heightDifference;

    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
