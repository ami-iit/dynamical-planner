/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <DynamicalPlannerPrivate/Costs/FootYawCost.h>
#include <iDynTree/Core/Utils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>

using namespace DynamicalPlanner::Private;

class FootYawCost::Implementation {
public:
    VariablesLabeller stateVariables;

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> desiredYaw;
    size_t bottomRightIndex, topRightIndex;
    double initialYaw;
    iDynTree::IndexRange bottomRightRange, topRightRange;

    iDynTree::VectorDynSize stateGradientBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};

FootYawCost::FootYawCost(const VariablesLabeller &stateVariables, const std::string &footName, const std::vector<iDynTree::Position> &pointsPosition)
    : iDynTree::optimalcontrol::Cost (footName + "FootYawCost")
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;

    m_pimpl->desiredYaw = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(0.0);

    assert(pointsPosition.size());

    iDynTree::Position centroid(0.0, 0.0, 0.0);

    for (auto& point: pointsPosition) {
        iDynTree::toEigen(centroid) += iDynTree::toEigen(point);
    }

    iDynTree::toEigen(centroid) /= pointsPosition.size();

    iDynTree::Vector3 relativePosition;
    double bottomRightValue = 0, topRightValue = 0;

    m_pimpl->bottomRightIndex = pointsPosition.size();
    m_pimpl->topRightIndex = pointsPosition.size();

    for (size_t i = 0; i < pointsPosition.size(); ++i) {
        relativePosition = pointsPosition[i] - centroid;
        if (relativePosition(1) < 0) {
            if (relativePosition(0) < 0) {
                if ((m_pimpl->bottomRightIndex == pointsPosition.size()) || (relativePosition(0) * relativePosition(1) > bottomRightValue)) {
                    m_pimpl->bottomRightIndex = i;
                    bottomRightValue = relativePosition(0) * relativePosition(1);
                }
            } else if (relativePosition(0) > 0) {
                if ((m_pimpl->topRightIndex == pointsPosition.size()) || (-relativePosition(0) * relativePosition(1) > topRightValue)) {
                    m_pimpl->topRightIndex = i;
                    topRightValue = -relativePosition(0) * relativePosition(1);
                }
            }
        }
    }
    assert(m_pimpl->bottomRightIndex < pointsPosition.size());
    assert(m_pimpl->topRightIndex < pointsPosition.size());
    assert(m_pimpl->bottomRightIndex != m_pimpl->topRightIndex);

    iDynTree::Vector3 initialDistance = pointsPosition[m_pimpl->topRightIndex] - pointsPosition[m_pimpl->bottomRightIndex];

    m_pimpl->initialYaw = std::atan2(initialDistance(1), initialDistance(0));

    m_pimpl->bottomRightRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(m_pimpl->bottomRightIndex));
    assert(m_pimpl->bottomRightRange.isValid());
    m_pimpl->topRightRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(m_pimpl->topRightIndex));
    assert(m_pimpl->topRightRange.isValid());

    m_pimpl->stateHessianSparsity.addDenseBlock(static_cast<size_t>(m_pimpl->topRightRange.offset),
                                                static_cast<size_t>(m_pimpl->topRightRange.offset), 2, 2);

    m_pimpl->stateHessianSparsity.addDenseBlock(static_cast<size_t>(m_pimpl->bottomRightRange.offset),
                                                static_cast<size_t>(m_pimpl->bottomRightRange.offset), 2, 2);

    m_pimpl->stateHessianSparsity.addDenseBlock(static_cast<size_t>(m_pimpl->topRightRange.offset),
                                                static_cast<size_t>(m_pimpl->bottomRightRange.offset), 2, 2);

    m_pimpl->stateHessianSparsity.addDenseBlock(static_cast<size_t>(m_pimpl->bottomRightRange.offset),
                                                static_cast<size_t>(m_pimpl->topRightRange.offset), 2, 2);

    m_pimpl->controlHessianSparsity.clear();
    m_pimpl->mixedHessianSparsity.clear();
}

FootYawCost::~FootYawCost()
{ }

void FootYawCost::setDesiredYawTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> desiredYaw)
{
    assert(desiredYaw);

    m_pimpl->desiredYaw = desiredYaw;
}

bool FootYawCost::costEvaluation(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, double &costValue)
{
    m_pimpl->stateVariables = state;

    bool isValid = false;
    double desiredYaw = m_pimpl->desiredYaw->get(time, isValid) - m_pimpl->initialYaw;
    assert(isValid);

    iDynTree::Vector3 currentDistance;
    iDynTree::toEigen(currentDistance) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->topRightRange)) - iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->bottomRightRange));

    double yawError = currentDistance(1) * std::cos(desiredYaw) - currentDistance(0) * std::sin(desiredYaw);

    costValue = 0.5 * yawError * yawError;

    return true;
}

bool FootYawCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    bool isValid = false;
    double desiredYaw = m_pimpl->desiredYaw->get(time, isValid) - m_pimpl->initialYaw;
    assert(isValid);

    iDynTree::Vector3 currentDistance;
    iDynTree::toEigen(currentDistance) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->topRightRange)) - iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->bottomRightRange));

    double sy = std::sin(desiredYaw);
    double cy = std::cos(desiredYaw);

    double yawError = currentDistance(1) * cy - currentDistance(0) * sy;

    partialDerivative(static_cast<unsigned int>(m_pimpl->topRightRange.offset)) = - yawError * sy;
    partialDerivative(static_cast<unsigned int>(m_pimpl->topRightRange.offset) + 1) = yawError * cy;

    partialDerivative(static_cast<unsigned int>(m_pimpl->bottomRightRange.offset)) = yawError * sy;
    partialDerivative(static_cast<unsigned int>(m_pimpl->bottomRightRange.offset) + 1) = - yawError * cy;

    return true;
}

bool FootYawCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &/*partialDerivative*/)
{
    return true;
}

bool FootYawCost::costSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    bool isValid = false;
    double desiredYaw = m_pimpl->desiredYaw->get(time, isValid) - m_pimpl->initialYaw;
    assert(isValid);

    iDynTree::Vector3 currentDistance;
    iDynTree::toEigen(currentDistance) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->topRightRange)) - iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->bottomRightRange));

    double sy = std::sin(desiredYaw);
    double cy = std::cos(desiredYaw);
    unsigned int trIndex = static_cast<unsigned int>(m_pimpl->topRightRange.offset);
    unsigned int brIndex = static_cast<unsigned int>(m_pimpl->bottomRightRange.offset);

    Eigen::Matrix<double, 2, 2> hessianBlock;
    hessianBlock(0,0) = sy * sy;
    hessianBlock(0,1) = -sy * cy;
    hessianBlock(1,0) = hessianBlock(0,1);
    hessianBlock(1,1) = cy * cy;

    iDynTree::toEigen(partialDerivative).block<2,2>(trIndex, trIndex) = hessianBlock;
    iDynTree::toEigen(partialDerivative).block<2,2>(brIndex, brIndex) = hessianBlock;
    iDynTree::toEigen(partialDerivative).block<2,2>(trIndex, brIndex) = -hessianBlock;
    iDynTree::toEigen(partialDerivative).block<2,2>(brIndex, trIndex) = -hessianBlock;


    return true;
}

bool FootYawCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FootYawCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FootYawCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool FootYawCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool FootYawCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
