/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/ForceMeanCost.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <vector>
#include <cassert>


using namespace DynamicalPlanner::Private;

class ForceMeanCost::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;

    iDynTree::IndexRange forcePointRange;
    std::vector<iDynTree::IndexRange> otherPointsRanges;

    iDynTree::Vector3 differenceFromMean, sumOfForces;
    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

};

ForceMeanCost::ForceMeanCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                             const std::string &footName, size_t contactIndex)
    : iDynTree::optimalcontrol::Cost ("ForceMean" + footName + std::to_string(contactIndex))
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    for (auto label : stateVariables.listOfLabels()) {
        if ((label.find(footName + "ForcePoint") != std::string::npos) && (label != (footName + "ForcePoint" + std::to_string(contactIndex)))) {
            m_pimpl->otherPointsRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->otherPointsRanges.back().isValid());
        }
    }

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();
}

bool ForceMeanCost::costEvaluation(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control,
                                   double &costValue)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sumOfForces = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    for (size_t i = i; i < m_pimpl->otherPointsRanges.size(); ++i) {
        iDynTree::toEigen(m_pimpl->sumOfForces) += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->otherPointsRanges[i]));
    }

    double numberOfPoints = m_pimpl->otherPointsRanges.size() + 1.0;

    iDynTree::toEigen(m_pimpl->differenceFromMean) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->forcePointRange)) -
            1.0/numberOfPoints * iDynTree::toEigen(m_pimpl->sumOfForces);

    costValue = 0.5 * iDynTree::toEigen(m_pimpl->differenceFromMean).squaredNorm();

    return true;
}

bool ForceMeanCost::costFirstPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

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
