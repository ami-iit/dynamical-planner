/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/Constraints/ContactForceControlConstraints.h>
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

    iDynTree::IndexRange positionPointRange, forcePointRange, forceControlRange;
    iDynTree::Vector3 pointPosition, pointForce, pointForceControl;

    iDynTree::VectorDynSize constraintValues;
    iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;

    iDynTree::optimalcontrol::SparsityStructure stateJacobianSparsity, controlJacobianSparsity;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;
};




ContactForceControlConstraints::ContactForceControlConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                               const std::string &footName, size_t contactIndex, const HyperbolicSecant &forceActivation,
                                                               double maximumNormalDerivative, double dissipationRatio)
    : iDynTree::optimalcontrol::Constraint (1, "ForceControlBounds" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_pimpl->footName = footName;
    m_pimpl->contactIndex = contactIndex;
    m_pimpl->activation = forceActivation;
    m_pimpl->maximumNormalDerivative = maximumNormalDerivative;
    m_pimpl->dissipationRatio = dissipationRatio;

    m_pimpl->positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
    assert(m_pimpl->positionPointRange.isValid());

    m_pimpl->forcePointRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->forceControlRange = controlVariables.getIndexRange(footName + "ForceControlPoint" + std::to_string(contactIndex));
    assert(m_pimpl->forcePointRange.isValid());

    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->constraintValues.resize(1);

    m_isLowerBounded = false;
    m_isUpperBounded = true;
    m_upperBound(0) = 0.0;

    m_pimpl->stateJacobianSparsity.clear();
    m_pimpl->controlJacobianSparsity.clear();

    size_t pzCol = static_cast<size_t>(m_pimpl->positionPointRange.offset + 2);
    size_t fzCol = static_cast<size_t>(m_pimpl->forcePointRange.offset + 2);
    size_t fzCtrlCol = static_cast<size_t>(m_pimpl->forceControlRange.offset + 2);

    m_pimpl->stateJacobianSparsity.add(0, pzCol);
    m_pimpl->stateJacobianSparsity.add(0, fzCol);

    m_pimpl->controlJacobianSparsity.add(0, fzCtrlCol);

    m_pimpl->stateHessianSparsity.add(pzCol, pzCol);
    m_pimpl->stateHessianSparsity.add(pzCol, fzCol);
    m_pimpl->stateHessianSparsity.add(fzCol, pzCol);

    m_pimpl->controlHessianSparsity.clear();
    m_pimpl->mixedHessianSparsity.clear();

}

ContactForceControlConstraints::~ContactForceControlConstraints()
{ }

bool ContactForceControlConstraints::evaluateConstraint(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);
    m_pimpl->pointForceControl = m_pimpl->controlVariables(m_pimpl->forceControlRange);
    double delta = m_pimpl->activation.eval(m_pimpl->pointPosition(2));
    double fz = m_pimpl->pointForce(2);
    double uz = m_pimpl->pointForceControl(2);

    m_pimpl->constraintValues(0) = uz - delta * m_pimpl->maximumNormalDerivative + (1- delta) * m_pimpl->dissipationRatio * fz;

    constraint = m_pimpl->constraintValues;

    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTState(double /*time*/, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
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

    m_pimpl->stateJacobianBuffer(0, pzCol) = -deltaDerivative * m_pimpl->maximumNormalDerivative -
        deltaDerivative * m_pimpl->dissipationRatio * fz;

    m_pimpl->stateJacobianBuffer(0, fzCol) = (1- delta) * m_pimpl->dissipationRatio;

    jacobian = m_pimpl->stateJacobianBuffer;
    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->controlJacobianBuffer(0, static_cast<unsigned int>(m_pimpl->forceControlRange.offset + 2)) = 1;

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
    stateSparsity = m_pimpl->stateJacobianSparsity;
    return true;
}

bool ContactForceControlConstraints::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlJacobianSparsity;
    return true;
}

bool ContactForceControlConstraints::constraintSecondPartialDerivativeWRTState(double /*time*/, const iDynTree::VectorDynSize &state,
                                                                               const iDynTree::VectorDynSize &/*control*/,
                                                                               const iDynTree::VectorDynSize &lambda,
                                                                               iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;

    m_pimpl->pointPosition = m_pimpl->stateVariables(m_pimpl->positionPointRange);
    m_pimpl->pointForce = m_pimpl->stateVariables(m_pimpl->forcePointRange);

    double fz = m_pimpl->pointForce(2);

    unsigned int pzIndex = static_cast<unsigned int>(m_pimpl->positionPointRange.offset + 2);
    unsigned int fzIndex = static_cast<unsigned int>(m_pimpl->forcePointRange.offset + 2);

    double maxDerivative = m_pimpl->maximumNormalDerivative;
    double deltaDerivative = m_pimpl->activation.evalDerivative(m_pimpl->pointPosition(2));
    double deltaDoubleDerivative = m_pimpl->activation.evalDoubleDerivative(m_pimpl->pointPosition(2));

    hessian(pzIndex, pzIndex) = -lambda(0) * (deltaDoubleDerivative * maxDerivative + deltaDoubleDerivative * m_pimpl->dissipationRatio * fz);
    hessian(pzIndex, fzIndex) = -lambda(0) * deltaDerivative * m_pimpl->dissipationRatio;
    hessian(fzIndex, pzIndex) = hessian(pzIndex, fzIndex);

    return true;
}

bool ContactForceControlConstraints::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                 const iDynTree::VectorDynSize &/*control*/,
                                                                                 const iDynTree::VectorDynSize &/*lambda*/,
                                                                                 iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ContactForceControlConstraints::constraintSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                      const iDynTree::VectorDynSize &/*control*/,
                                                                                      const iDynTree::VectorDynSize &/*lambda*/,
                                                                                      iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ContactForceControlConstraints::constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ContactForceControlConstraints::constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ContactForceControlConstraints::constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
