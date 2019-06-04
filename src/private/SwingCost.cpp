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
    iDynTree::Vector3 weights;

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};

SwingCost::SwingCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, const std::string &footName,
                     size_t contactIndex, double desiredSwingHeight, const iDynTree::Vector3 &weights)
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

    m_pimpl->weights = weights;

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    unsigned int uIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);

    m_pimpl->stateHessianSparsity.add(pzIndex, pzIndex);
    m_pimpl->controlHessianSparsity.addIdentityBlock(uIndex, uIndex, 3);
    m_pimpl->mixedHessianSparsity.addDenseBlock(pzIndex, uIndex, 1, 3);
}

SwingCost::~SwingCost()
{ }

bool SwingCost::costEvaluation(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    double zError = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    double costX = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) * zError;
    double costY = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) * zError;
    double costZ = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(2) * zError;

    costValue = 0.5 * (m_pimpl->weights(0) * costX * costX +
                       m_pimpl->weights(1) * costY * costY +
                       m_pimpl->weights(2) * costZ * costZ);
    return true;
}

bool SwingCost::costFirstPartialDerivativeWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    double zError = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    double costX = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) * zError;
    double costY = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) * zError;
    double costZ = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(2) * zError;


    m_pimpl->stateGradientBuffer(pzIndex) = m_pimpl->weights(0) * costX * m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) +
        m_pimpl->weights(1) * costY * m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) +
        m_pimpl->weights(2) * costZ * m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(2);

    partialDerivative = m_pimpl->stateGradientBuffer;

    return true;
}

bool SwingCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);
    unsigned int uzIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 2);


    double zError = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    double costX = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0) * zError;
    double costY = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1) * zError;
    double costZ = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(2) * zError;

    m_pimpl->controlGradientBuffer(uxIndex) = m_pimpl->weights(0) * costX * zError;
    m_pimpl->controlGradientBuffer(uyIndex) = m_pimpl->weights(1) * costY * zError;
    m_pimpl->controlGradientBuffer(uzIndex) = m_pimpl->weights(2) * costZ * zError;

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
    double uz = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(2);

    partialDerivative(pzIndex, pzIndex) = m_pimpl->weights(0) * ux * ux +
        m_pimpl->weights(1) * uy * uy +
        m_pimpl->weights(2) * uz * uz;

    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);
    unsigned int uzIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 2);
    double heightDifference = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);
    double squaredHeightDifference = heightDifference * heightDifference;

    partialDerivative(uxIndex, uxIndex) = m_pimpl->weights(0) * squaredHeightDifference;
    partialDerivative(uyIndex, uyIndex) = m_pimpl->weights(1) * squaredHeightDifference;
    partialDerivative(uzIndex, uzIndex) = m_pimpl->weights(2) * squaredHeightDifference;

    return true;
}

bool SwingCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->pointPositionRange.offset + 2);
    unsigned int uxIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset);
    unsigned int uyIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 1);
    unsigned int uzIndex = static_cast<unsigned int>(m_pimpl->pointVelocityRange.offset + 2);

    double ux = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(0);
    double uy = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(1);
    double uz = m_pimpl->controlVariables(m_pimpl->pointVelocityRange)(2);

    double heightDifference = (m_pimpl->stateVariables(m_pimpl->pointPositionRange)(2) - m_pimpl->swingHeight);

    partialDerivative(pzIndex, uxIndex) = 2.0 * m_pimpl->weights(0) * ux * heightDifference;
    partialDerivative(pzIndex, uyIndex) = 2.0 * m_pimpl->weights(1) * uy * heightDifference;
    partialDerivative(pzIndex, uzIndex) = 2.0 * m_pimpl->weights(2) * uz * heightDifference;


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
