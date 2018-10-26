/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/ContactForceControlConstraints.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/MatrixFixSize.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class ContactForceControlConstraints::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    std::string footName;
    size_t contactIndex;
    HyperbolicSecant activation;
    double maximumNormalDerivative;
    double dissipationRatio;
    double deactivationTime;

    iDynTree::IndexRange positionPointRange, forcePointRange, forceControlRange;
    iDynTree::Vector3 pointPosition, pointForce, pointForceControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;
};




ContactForceControlConstraints::ContactForceControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                               const std::string &footName, size_t contactIndex, const HyperbolicSecant &forceActivation,
                                                               double maximumNormalDerivative, double dissipationRatio, double deactivationTime)
    : iDynTree::optimalcontrol::Constraint (2, "ForceControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->activation = forceActivation;
    m_pimpl->maximumNormalDerivative = maximumNormalDerivative;
    m_pimpl->dissipationRatio = dissipationRatio;
    m_pimpl->deactivationTime = deactivationTime;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->forceControlRange = controlVariables.getIndexRange(footName + "ForceControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(2, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(2, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(2);

    m_isLowerBounded = true;
    m_isUpperBounded = false;
    m_lowerBound.zero();

    m_pimpl->stateSparsity.clear();
    m_pimpl->controlSparsity.clear();

    size_t pzCol = static_cast<size_t>(m_pimpl->positionPointRange.offset + 2);
    size_t fzCol = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    size_t fzCtrlCol = static_cast<size_t>(m_pimpl->forceControlRange.offset + 2);

    m_pimpl->stateSparsity.addDenseBlock(0, pzCol, 2, 1);
    m_pimpl->stateSparsity.addDenseBlock(0, fzCol, 2, 1);

    m_pimpl->controlSparsity.addDenseBlock(0, fzCtrlCol, 2, 1);


}

ContactForceControlConstraints::~ContactForceControlConstraints()
{ }

bool ContactForceControlConstraints::evaluateConstraint(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);
    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double fz = m_pimpl->pointForce(2);
    double uz = m_pimpl->pointForceControl(2);

    if (time > m_pimpl->deactivationTime) {
        m_pimpl->constraintValues(0) = - (1- delta) * m_pimpl->dissipationRatio * fz - uz;

        m_pimpl->constraintValues(1) = uz + (1- delta) * m_pimpl->dissipationRatio * fz;
    } else {
        m_pimpl->constraintValues(0) = delta * m_pimpl->maximumNormalDerivative - (1- delta) * m_pimpl->dissipationRatio * fz - uz;

        m_pimpl->constraintValues(1) = uz + delta * m_pimpl->maximumNormalDerivative + (1- delta) * m_pimpl->dissipationRatio * fz;
    }

    constraint = m_pimpl->constraintValues;

    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);

    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double deltaDerivative = m_pimpl->activation.evalDerivative(m_pimpl->pointPosition(2));
    double fz = m_pimpl->pointForce(2);
    unsigned int pzCol = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);
    unsigned int fzCol = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    if (time > m_pimpl->deactivationTime) {
        m_pimpl->stateJacobianBuffer(0, pzCol) = +deltaDerivative * m_pimpl->dissipationRatio * fz;

        m_pimpl->stateJacobianBuffer(0, fzCol) = -(1- delta) * m_pimpl->dissipationRatio;

        m_pimpl->stateJacobianBuffer(1, pzCol) = -deltaDerivative * m_pimpl->dissipationRatio * fz;

        m_pimpl->stateJacobianBuffer(1, fzCol) = (1- delta) * m_pimpl->dissipationRatio;
    } else {
        m_pimpl->stateJacobianBuffer(0, pzCol) = deltaDerivative * m_pimpl->maximumNormalDerivative +
                deltaDerivative * m_pimpl->dissipationRatio * fz;

        m_pimpl->stateJacobianBuffer(0, fzCol) = -(1- delta) * m_pimpl->dissipationRatio;

        m_pimpl->stateJacobianBuffer(1, pzCol) = deltaDerivative * m_pimpl->maximumNormalDerivative -
                deltaDerivative * m_pimpl->dissipationRatio * fz;

        m_pimpl->stateJacobianBuffer(1, fzCol) = (1- delta) * m_pimpl->dissipationRatio;
    }

    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->controlJacobianBuffer(0, static_cast<unsigned int>(m_pimpl->forceControlRange.offset + 2)) = -1;

    m_pimpl->controlJacobianBuffer(1, static_cast<unsigned int>(m_pimpl->forceControlRange.offset + 2)) = 1;

    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ContactForceControlConstraints::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ContactForceControlConstraints::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool ContactForceControlConstraints::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}
