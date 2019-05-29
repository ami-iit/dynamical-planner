/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/FeetDistanceCost.h>
#include <iDynTree/Core/Utils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <vector>
#include <cassert>
#include <iostream>

using namespace DynamicalPlanner::Private;

class FeetDistanceCost::Implementation {
public:
    VariablesLabeller stateVariables;

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> desiredDistance;

    std::vector<iDynTree::IndexRange> leftPointRanges, rightPointRanges;
    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

    iDynTree::MatrixDynSize stateHessian;
};

FeetDistanceCost::FeetDistanceCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables)
    :iDynTree::optimalcontrol::Cost ("FeetDistance")
      , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;

    m_pimpl->desiredDistance = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(0.0);

    std::vector<std::pair<iDynTree::IndexRange, bool>> allPoints;

    for (auto& label : stateVariables.listOfLabels()) {
        if (label.find("LeftPositionPoint") != std::string::npos) {
            m_pimpl->leftPointRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->leftPointRanges.back().isValid());
            allPoints.push_back(std::make_pair(stateVariables.getIndexRange(label), true));
        }
        if (label.find("RightPositionPoint") != std::string::npos) {
            m_pimpl->rightPointRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->rightPointRanges.back().isValid());
            allPoints.push_back(std::make_pair(stateVariables.getIndexRange(label), false));
        }
    }
    assert(m_pimpl->leftPointRanges.size() && m_pimpl->rightPointRanges.size());

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();

    m_pimpl->stateHessian.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateHessian.zero();

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(m_pimpl->stateHessian);
    double nLeft = m_pimpl->leftPointRanges.size();
    double nRight = m_pimpl->rightPointRanges.size();

    auto diagonalValue = [nLeft, nRight, allPoints](size_t i, size_t j ) {
        bool iIsLeft = allPoints[i].second;
        bool jIsLeft = allPoints[j].second;

        if (iIsLeft) {
            if (jIsLeft) {
                return nLeft * nLeft;
            } else {
                return -nLeft * nRight;
            }
        } else {
            if (jIsLeft) {
                return -nLeft * nRight;
            } else {
                return nRight * nRight;
            }
        }
    };

    for (size_t i = 0; i < allPoints.size(); ++i) {
        for (size_t j = i; j < allPoints.size(); ++j) {

            unsigned int iOffset = static_cast<unsigned int>(allPoints[i].first.offset);
            unsigned int jOffset = static_cast<unsigned int>(allPoints[j].first.offset);

            m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(allPoints[i].first.offset),
                                                           static_cast<size_t>(allPoints[j].first.offset), 3);

            hessianMap.block<3,3>(iOffset, jOffset).setIdentity();
            hessianMap.block<3,3>(iOffset, jOffset) /= diagonalValue(i, j);

            if (i != j) {
                m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(allPoints[j].first.offset),
                                                               static_cast<size_t>(allPoints[i].first.offset), 3);
                hessianMap.block<3,3>(jOffset, iOffset) = hessianMap.block<3,3>(iOffset, jOffset);
            }
        }
    }

    m_pimpl->mixedHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();

}

FeetDistanceCost::~FeetDistanceCost()
{

}

bool FeetDistanceCost::costEvaluation(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, double &costValue)
{
    m_pimpl->stateVariables = state;

    Eigen::Vector3d leftPosition, rightPosition, distance;

    leftPosition = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->leftPointRanges.front()));
    for (size_t i = 1; i < m_pimpl->leftPointRanges.size(); ++i) {
        leftPosition += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->leftPointRanges[i]));
    }

    rightPosition = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->rightPointRanges.front()));
    for (size_t i = 1; i < m_pimpl->rightPointRanges.size(); ++i) {
        rightPosition += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->rightPointRanges[i]));
    }

    distance = leftPosition/m_pimpl->leftPointRanges.size() - rightPosition/m_pimpl->rightPointRanges.size();

    costValue = 0.5 * distance.squaredNorm();

    return true;
}

bool FeetDistanceCost::costFirstPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    Eigen::Vector3d leftPosition, rightPosition, distance;

    leftPosition = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->leftPointRanges.front()));
    for (size_t i = 1; i < m_pimpl->leftPointRanges.size(); ++i) {
        leftPosition += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->leftPointRanges[i]));
    }

    rightPosition = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->rightPointRanges.front()));
    for (size_t i = 1; i < m_pimpl->rightPointRanges.size(); ++i) {
        rightPosition += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->rightPointRanges[i]));
    }

    distance = leftPosition/m_pimpl->leftPointRanges.size() - rightPosition/m_pimpl->rightPointRanges.size();

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(partialDerivative);

    double nLeft = m_pimpl->leftPointRanges.size();

    for (auto& point : m_pimpl->leftPointRanges) {
        gradientMap.segment<3>(point.offset) = distance/nLeft;
    }

    double nRight = m_pimpl->rightPointRanges.size();

    for (auto& point : m_pimpl->rightPointRanges) {
        gradientMap.segment<3>(point.offset) = -distance/nRight;
    }

    return true;
}

bool FeetDistanceCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &/*partialDerivative*/)
{
    return true;
}

bool FeetDistanceCost::costSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->stateHessian;

    return true;
}

bool FeetDistanceCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                  const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FeetDistanceCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                       const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool FeetDistanceCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool FeetDistanceCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool FeetDistanceCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
