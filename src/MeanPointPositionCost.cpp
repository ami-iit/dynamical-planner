/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/MeanPointPositionCost.h>
#include <iDynTree/Core/Utils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <vector>
#include <cassert>
#include <iostream>

using namespace DynamicalPlanner::Private;

class MeanPointPositionCost::Implementation {
public:
    VariablesLabeller stateVariables;

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> desiredPosition;

    std::vector<iDynTree::IndexRange> pointRanges;
    iDynTree::Vector3 distanceFromTarget;
    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    iDynTree::optimalcontrol::TimeRange horizon;
    double increaseFactor;
};


MeanPointPositionCost::MeanPointPositionCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables)
    :iDynTree::optimalcontrol::Cost ("MeanPointPosition")
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;

    m_pimpl->desiredPosition = std::make_shared<iDynTree::optimalcontrol::TimeInvariantPosition>(iDynTree::Position::Zero());

    for (auto& label : stateVariables.listOfLabels()) {
        if (label.find("PositionPoint") != std::string::npos) {
            m_pimpl->pointRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->pointRanges.back().isValid());
        }
    }
    assert(m_pimpl->pointRanges.size());

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();

    m_pimpl->increaseFactor = 0.0;
    m_pimpl->horizon.setTimeInterval(0.0, 1.0);

}

MeanPointPositionCost::~MeanPointPositionCost()
{

}

bool MeanPointPositionCost::setDesiredPositionTrajectory(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> desiredPosition)
{
    if (!desiredPosition) {
        return false;
    }
    m_pimpl->desiredPosition = desiredPosition;
    return true;
}

bool MeanPointPositionCost::setTimePenalty(const iDynTree::optimalcontrol::TimeRange &horizon, double increaseFactor)
{
    if (!horizon.isValid() || horizon.isInstant()) {
        std::cerr << "[ERROR][MeanPointPositionCost::setTimePenalty] The horizon is " << (horizon.isInstant() ? "a single instant." : "not valid.") << std::endl;
        return false;
    }

    if (increaseFactor < 0) {
        std::cerr << "[ERROR][MeanPointPositionCost::setTimePenalty] The increase factor is supposed to be non-negative." << std::endl;
        return false;
    }

    m_pimpl->horizon = horizon;
    m_pimpl->increaseFactor = increaseFactor;

    return true;
}

bool MeanPointPositionCost::costEvaluation(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &,
                                           double &costValue)
{
    m_pimpl->stateVariables = state;

    m_pimpl->distanceFromTarget = m_pimpl->stateVariables(m_pimpl->pointRanges[0]);
    for (size_t i = 1; i < m_pimpl->pointRanges.size(); ++i) {
        iDynTree::toEigen(m_pimpl->distanceFromTarget) += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->pointRanges[i]));
    }
    iDynTree::toEigen(m_pimpl->distanceFromTarget) /= m_pimpl->pointRanges.size();

    bool isValid;
    const iDynTree::Position& desiredPosition = m_pimpl->desiredPosition->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][MeanPointPositionCost::costEvaluation] Unable to retrieve a valid position at time " << time
                  << "." << std::endl;
        return false;
    }

    iDynTree::toEigen(m_pimpl->distanceFromTarget) -= iDynTree::toEigen(desiredPosition);

    double timeRatio = (time - m_pimpl->horizon.initTime()) / (m_pimpl->horizon.endTime() - m_pimpl->horizon.initTime());
    double timeWeight = m_pimpl->increaseFactor * std::fabs(timeRatio) + 1.0;
    double timeWeightSquared = timeWeight * timeWeight;

    costValue = 0.5 * timeWeightSquared * iDynTree::toEigen(m_pimpl->distanceFromTarget).squaredNorm();
    return true;
}

bool MeanPointPositionCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                               const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    double numberOfPointsInverse = 1/static_cast<double>(m_pimpl->pointRanges.size());

    m_pimpl->distanceFromTarget = m_pimpl->stateVariables(m_pimpl->pointRanges[0]);
    for (size_t i = 1; i < m_pimpl->pointRanges.size(); ++i) {
        iDynTree::toEigen(m_pimpl->distanceFromTarget) += iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->pointRanges[i]));
    }
    iDynTree::toEigen(m_pimpl->distanceFromTarget) *= numberOfPointsInverse;

    bool isValid;
    const iDynTree::Position& desiredPosition = m_pimpl->desiredPosition->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][MeanPointPositionCost::costFirstPartialDerivativeWRTState] Unable to retrieve a valid position at time " << time
                  << "." << std::endl;
        return false;
    }

    iDynTree::toEigen(m_pimpl->distanceFromTarget) -= iDynTree::toEigen(desiredPosition);

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);

    double timeRatio = (time - m_pimpl->horizon.initTime()) / (m_pimpl->horizon.endTime() - m_pimpl->horizon.initTime());
    double timeWeight = m_pimpl->increaseFactor * std::fabs(timeRatio) + 1.0;
    double timeWeightSquared = timeWeight * timeWeight;

    for (auto& point : m_pimpl->pointRanges) {
        gradientMap.segment<3>(point.offset) = timeWeightSquared *  numberOfPointsInverse * iDynTree::toEigen(m_pimpl->distanceFromTarget);
    }

    partialDerivative = m_pimpl->stateGradientBuffer;
    return true;
}

bool MeanPointPositionCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &,
                                                                 iDynTree::VectorDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->controlGradientBuffer;
    return true;
}
