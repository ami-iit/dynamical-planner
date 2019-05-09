/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/ForceMeanCost.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <vector>
#include <cassert>


using namespace DynamicalPlanner::Private;

class ForceMeanCost::Implementation {
public:
    VariablesLabeller stateVariables;

    std::string footName;
    size_t contactIndex;

    iDynTree::IndexRange forcePointRange;
    std::vector<iDynTree::IndexRange> otherPointsRanges;

    iDynTree::Vector3 differenceFromMean, sumOfForces;
    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

};

ForceMeanCost::ForceMeanCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                             const std::string &footName, size_t contactIndex)
    : iDynTree::optimalcontrol::Cost ("ForceMean" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    for (auto& label : stateVariables.listOfLabels()) {
        if ((label.find(footName + "ForcePoint") != std::string::npos) && (label != (footName + "ForcePoint" + std::to_string(contactIndex)))) {
            m_pimpl->otherPointsRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->otherPointsRanges.back().isValid());
        }
    }

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();

    m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->forcePointRange.offset),
                                                   static_cast<size_t>(m_pimpl->forcePointRange.offset), 3);

    for (auto& force : m_pimpl->otherPointsRanges) {
        m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->forcePointRange.offset),
                                                       static_cast<size_t>(force.offset), 3);
        m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(force.offset),
                                                       static_cast<size_t>(m_pimpl->forcePointRange.offset), 3);
    }

    for (size_t i = 0; i < m_pimpl->otherPointsRanges.size(); ++i) {
        for (size_t j = i; j < m_pimpl->otherPointsRanges.size(); ++j) {

            m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->otherPointsRanges[i].offset),
                                                           static_cast<size_t>(m_pimpl->otherPointsRanges[j].offset), 3);

            if (i != j) {
                m_pimpl->stateHessianSparsity.addIdentityBlock(static_cast<size_t>(m_pimpl->otherPointsRanges[j].offset),
                                                               static_cast<size_t>(m_pimpl->otherPointsRanges[i].offset), 3);
            }
        }
    }

    m_pimpl->mixedHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();
}

ForceMeanCost::~ForceMeanCost()
{ }

bool ForceMeanCost::costEvaluation(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &,
                                   double &costValue)
{
    m_pimpl->stateVariables = state;

    m_pimpl->sumOfForces = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    for (auto force : m_pimpl->otherPointsRanges) {
        iDynTree::toEigen(m_pimpl->sumOfForces) += iDynTree::toEigen(m_pimpl->stateVariables(force));
    }

    double numberOfPoints = m_pimpl->otherPointsRanges.size() + 1.0;

    iDynTree::toEigen(m_pimpl->differenceFromMean) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->forcePointRange)) -
            1.0/numberOfPoints * iDynTree::toEigen(m_pimpl->sumOfForces);

    costValue = 0.5 * iDynTree::toEigen(m_pimpl->differenceFromMean).squaredNorm();

    return true;
}

bool ForceMeanCost::costFirstPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    m_pimpl->sumOfForces = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    for (auto force : m_pimpl->otherPointsRanges) {
        iDynTree::toEigen(m_pimpl->sumOfForces) += iDynTree::toEigen(m_pimpl->stateVariables(force));
    }

    double numberOfPoints = m_pimpl->otherPointsRanges.size() + 1.0;
    double numberOfPointsInverse = 1.0 / numberOfPoints;

    iDynTree::toEigen(m_pimpl->differenceFromMean) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->forcePointRange)) -
            1.0/numberOfPoints * iDynTree::toEigen(m_pimpl->sumOfForces);

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);

    gradientMap.segment<3>(m_pimpl->forcePointRange.offset) = (1.0 - numberOfPointsInverse) * iDynTree::toEigen(m_pimpl->differenceFromMean);

    for (auto force : m_pimpl->otherPointsRanges) {
        gradientMap.segment<3>(force.offset) = -numberOfPointsInverse * iDynTree::toEigen(m_pimpl->differenceFromMean);
    }

    partialDerivative = m_pimpl->stateGradientBuffer;

    return true;

}

bool ForceMeanCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->controlGradientBuffer;
    return true;
}

bool ForceMeanCost::costSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                        const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    double numberOfPoints = m_pimpl->otherPointsRanges.size() + 1.0;
    double numberOfPointsInverse = 1.0 / numberOfPoints;

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(partialDerivative);

    hessianMap.block<3,3>(m_pimpl->forcePointRange.offset, m_pimpl->forcePointRange.offset).setIdentity();
    hessianMap.block<3,3>(m_pimpl->forcePointRange.offset, m_pimpl->forcePointRange.offset) *= (1.0 - numberOfPointsInverse) *
        (1.0 - numberOfPointsInverse);

    for (auto& force : m_pimpl->otherPointsRanges) {
        hessianMap.block<3,3>(m_pimpl->forcePointRange.offset, force.offset).setIdentity();
        hessianMap.block<3,3>(m_pimpl->forcePointRange.offset, force.offset) *= -numberOfPointsInverse * (1.0 - numberOfPointsInverse);
        hessianMap.block<3,3>(force.offset, m_pimpl->forcePointRange.offset) =
            hessianMap.block<3,3>(m_pimpl->forcePointRange.offset, force.offset);
    }

    for (size_t i = 0; i < m_pimpl->otherPointsRanges.size(); ++i) {
        for (size_t j = i; j < m_pimpl->otherPointsRanges.size(); ++j) {
            hessianMap.block<3,3>(m_pimpl->otherPointsRanges[i].offset, m_pimpl->otherPointsRanges[j].offset).setIdentity();
            hessianMap.block<3,3>(m_pimpl->otherPointsRanges[i].offset, m_pimpl->otherPointsRanges[j].offset) *= numberOfPointsInverse *
                numberOfPointsInverse;

            if (i != j) {
                hessianMap.block<3,3>(m_pimpl->otherPointsRanges[j].offset, m_pimpl->otherPointsRanges[i].offset) =
                    hessianMap.block<3,3>(m_pimpl->otherPointsRanges[i].offset, m_pimpl->otherPointsRanges[j].offset);
            }
        }
    }

    return true;
}

bool ForceMeanCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                          const iDynTree::VectorDynSize &/*control*/,
                                                          iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool ForceMeanCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                               const iDynTree::VectorDynSize &/*control*/,
                                                               iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool ForceMeanCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ForceMeanCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ForceMeanCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
