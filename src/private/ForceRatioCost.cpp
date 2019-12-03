/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Costs/ForceRatioCost.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <vector>
#include <cassert>


using namespace DynamicalPlanner::Private;

class ForceRatioCost::Implementation {
public:
    VariablesLabeller stateVariables;

    std::string footName;
    size_t contactIndex;

    iDynTree::IndexRange forcePointRange;
    std::vector<iDynTree::IndexRange> otherPointsRanges;

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;

    double desiredRatio;

    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

};

ForceRatioCost::ForceRatioCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                 const std::string &footName, size_t contactIndex)
    : iDynTree::optimalcontrol::Cost ("ForceRatio" + footName + std::to_string(contactIndex))
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

    m_pimpl->stateHessianSparsity.add(static_cast<size_t>(m_pimpl->forcePointRange.offset) + 2,
                                      static_cast<size_t>(m_pimpl->forcePointRange.offset) + 2);

    for (auto& force : m_pimpl->otherPointsRanges) {
        m_pimpl->stateHessianSparsity.add(static_cast<size_t>(m_pimpl->forcePointRange.offset) + 2,
                                          static_cast<size_t>(force.offset) + 2);
        m_pimpl->stateHessianSparsity.add(static_cast<size_t>(force.offset) + 2,
                                          static_cast<size_t>(m_pimpl->forcePointRange.offset) + 2);
    }

    for (size_t i = 0; i < m_pimpl->otherPointsRanges.size(); ++i) {
        for (size_t j = i; j < m_pimpl->otherPointsRanges.size(); ++j) {

            m_pimpl->stateHessianSparsity.add(static_cast<size_t>(m_pimpl->otherPointsRanges[i].offset) + 2,
                                              static_cast<size_t>(m_pimpl->otherPointsRanges[j].offset) + 2);

            if (i != j) {
                m_pimpl->stateHessianSparsity.add(static_cast<size_t>(m_pimpl->otherPointsRanges[j].offset) + 2,
                                                  static_cast<size_t>(m_pimpl->otherPointsRanges[i].offset) + 2);
            }
        }
    }

    m_pimpl->desiredRatio = 1.0 / (m_pimpl->otherPointsRanges.size() + 1.0);

    m_pimpl->mixedHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();
}

ForceRatioCost::~ForceRatioCost()
{ }

void ForceRatioCost::setDesiredRatio(double ratio)
{
    m_pimpl->desiredRatio = ratio;
}

bool ForceRatioCost::costEvaluation(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &,
                                     double &costValue)
{
    m_pimpl->stateVariables = state;

    double sumOfForces = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);
    for (auto force : m_pimpl->otherPointsRanges) {
        sumOfForces += m_pimpl->stateVariables(force)(2);
    }

    double differenceFromRatio = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2) - m_pimpl->desiredRatio * sumOfForces;

    costValue = 0.5 * differenceFromRatio * differenceFromRatio;

    return true;
}

bool ForceRatioCost::costFirstPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    m_pimpl->stateVariables = state;

    double sumOfForces = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2);
    for (auto force : m_pimpl->otherPointsRanges) {
        sumOfForces += m_pimpl->stateVariables(force)(2);
    }

    double differenceFromRatio = m_pimpl->stateVariables(m_pimpl->forcePointRange)(2) - m_pimpl->desiredRatio * sumOfForces;

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);

    gradientMap(m_pimpl->forcePointRange.offset + 2) = (1.0 - m_pimpl->desiredRatio) * differenceFromRatio;

    for (auto force : m_pimpl->otherPointsRanges) {
        gradientMap(force.offset + 2) = -m_pimpl->desiredRatio * differenceFromRatio;
    }

    partialDerivative = m_pimpl->stateGradientBuffer;

    return true;

}

bool ForceRatioCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->controlGradientBuffer;
    return true;
}

bool ForceRatioCost::costSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                        const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &partialDerivative)
{
    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(partialDerivative);

    hessianMap(m_pimpl->forcePointRange.offset + 2, m_pimpl->forcePointRange.offset + 2) = (1.0 - m_pimpl->desiredRatio) * (1.0 - m_pimpl->desiredRatio);

    for (auto& force : m_pimpl->otherPointsRanges) {
        hessianMap(m_pimpl->forcePointRange.offset + 2, force.offset + 2) = - m_pimpl->desiredRatio * (1.0 - m_pimpl->desiredRatio);
        hessianMap(force.offset + 2, m_pimpl->forcePointRange.offset + 2) = hessianMap(m_pimpl->forcePointRange.offset + 2, force.offset + 2);
    }

    for (size_t i = 0; i < m_pimpl->otherPointsRanges.size(); ++i) {
        for (size_t j = i; j < m_pimpl->otherPointsRanges.size(); ++j) {
            hessianMap(m_pimpl->otherPointsRanges[i].offset + 2, m_pimpl->otherPointsRanges[j].offset + 2) = m_pimpl->desiredRatio * m_pimpl->desiredRatio;

            if (i != j) {
                hessianMap(m_pimpl->otherPointsRanges[j].offset + 2, m_pimpl->otherPointsRanges[i].offset + 2) = hessianMap(m_pimpl->otherPointsRanges[i].offset + 2, m_pimpl->otherPointsRanges[j].offset + 2);
            }
        }
    }

    return true;
}

bool ForceRatioCost::costSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                          const iDynTree::VectorDynSize &/*control*/,
                                                          iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool ForceRatioCost::costSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                               const iDynTree::VectorDynSize &/*control*/,
                                                               iDynTree::MatrixDynSize &/*partialDerivative*/)
{
    return true;
}

bool ForceRatioCost::costSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ForceRatioCost::costSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ForceRatioCost::costSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
