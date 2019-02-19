/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/MeanPointPositionCost.h>
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

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> timeVaryingWeight;
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

    m_pimpl->timeVaryingWeight = std::make_shared<iDynTree::optimalcontrol::TimeInvariantDouble>(1.0);

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

void MeanPointPositionCost::setTimeVaryingWeight(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> timeVaryingWeight)
{
    assert(timeVaryingWeight);

    m_pimpl->timeVaryingWeight = timeVaryingWeight;
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

    const double& timeWeight = m_pimpl->timeVaryingWeight->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][MeanPointPositionCost::costEvaluation] Unable to retrieve a valid timeVaryingWeight at time " << time
                  << "." << std::endl;
        return false;
    }

    iDynTree::toEigen(m_pimpl->distanceFromTarget) -= iDynTree::toEigen(desiredPosition);

    costValue = 0.5 * timeWeight * iDynTree::toEigen(m_pimpl->distanceFromTarget).squaredNorm();
    return true;
}

bool MeanPointPositionCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                               const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    double numberOfPointsInverse = 1.0/static_cast<double>(m_pimpl->pointRanges.size());

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

    const double& timeWeight = m_pimpl->timeVaryingWeight->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][MeanPointPositionCost::costEvaluation] Unable to retrieve a valid timeVaryingWeight at time " << time
                  << "." << std::endl;
        return false;
    }

    iDynTree::toEigen(m_pimpl->distanceFromTarget) -= iDynTree::toEigen(desiredPosition);

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);

    for (auto& point : m_pimpl->pointRanges) {
        gradientMap.segment<3>(point.offset) = timeWeight *  numberOfPointsInverse * iDynTree::toEigen(m_pimpl->distanceFromTarget);
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

bool MeanPointPositionCost::costSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    double numberOfPointsInverse = 1.0/static_cast<double>(m_pimpl->pointRanges.size());

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(partialDerivative);

    bool isValid;
    const double& timeWeight = m_pimpl->timeVaryingWeight->get(time, isValid);

    if (!isValid) {
        std::cerr << "[ERROR][MeanPointPositionCost::costEvaluation] Unable to retrieve a valid timeVaryingWeight at time " << time
                  << "." << std::endl;
        return false;
    }


    for (size_t i = 0; i < m_pimpl->pointRanges.size(); ++i) {
        for (size_t j = i; j < m_pimpl->pointRanges.size(); ++j) {
            hessianMap.block<3,3>(m_pimpl->pointRanges[i].offset, m_pimpl->pointRanges[j].offset).setIdentity();
            hessianMap.block<3,3>(m_pimpl->pointRanges[i].offset, m_pimpl->pointRanges[j].offset) *= timeWeight * numberOfPointsInverse * numberOfPointsInverse;

            if (i != j) {
                hessianMap.block<3,3>(m_pimpl->pointRanges[j].offset, m_pimpl->pointRanges[i].offset) =
                    hessianMap.block<3,3>(m_pimpl->pointRanges[i].offset, m_pimpl->pointRanges[j].offset);
            }
        }
    }

    return true;
}

bool MeanPointPositionCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                  const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool MeanPointPositionCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                       const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}
