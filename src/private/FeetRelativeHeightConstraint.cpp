/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <DynamicalPlannerPrivate/Constraints/FeetRelativeHeightConstraint.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class FeetRelativeHeightConstraint::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    iDynTree::MatrixDynSize relativeJacobianBuffer, stateJacobianBuffer, controlJacobianBuffer;

    std::vector<iDynTree::IndexRange> leftPointRanges, rightPointRanges;


    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

};


FeetRelativeHeightConstraint::FeetRelativeHeightConstraint(const VariablesLabeller& stateVariables,
                                                           const VariablesLabeller& controlVariables,
                                                           double maxRelativeHeight)
    : iDynTree::optimalcontrol::Constraint (1, "MaxRelativeHeight")
      , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    for (auto& label : stateVariables.listOfLabels()) {
        if (label.find("LeftPositionPoint") != std::string::npos) {
            m_pimpl->leftPointRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->leftPointRanges.back().isValid());
            m_pimpl->stateSparsity.add(0, static_cast<size_t>(m_pimpl->leftPointRanges.back().offset)+2);
        }
    }
    assert(m_pimpl->leftPointRanges.size());

    for (auto& label : stateVariables.listOfLabels()) {
        if (label.find("RightPositionPoint") != std::string::npos) {
            m_pimpl->rightPointRanges.push_back(stateVariables.getIndexRange(label));
            assert(m_pimpl->rightPointRanges.back().isValid());
            m_pimpl->stateSparsity.add(0, static_cast<size_t>(m_pimpl->rightPointRanges.back().offset)+2);
        }
    }
    assert(m_pimpl->rightPointRanges.size());

    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    for (auto& range : m_pimpl->leftPointRanges) {
        m_pimpl->stateJacobianBuffer(0, static_cast<unsigned int>(range.offset) + 2) = 1.0/m_pimpl->leftPointRanges.size();
    }

    for (auto& range : m_pimpl->rightPointRanges) {
        m_pimpl->stateJacobianBuffer(0, static_cast<unsigned int>(range.offset) + 2) = -1.0/m_pimpl->rightPointRanges.size();
    }

    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();


    m_lowerBound(0) = -maxRelativeHeight;
    m_upperBound(0) = maxRelativeHeight;

    m_isLowerBounded = true;
    m_isUpperBounded = true;


    m_pimpl->controlSparsity.clear();
    m_pimpl->stateHessianSparsity.clear();
    m_pimpl->controlHessianSparsity.clear();
    m_pimpl->mixedHessianSparsity.clear();

}

FeetRelativeHeightConstraint::~FeetRelativeHeightConstraint()
{ }

bool FeetRelativeHeightConstraint::evaluateConstraint(double , const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;

    double leftHeight = 0.0;
    for (auto& range : m_pimpl->leftPointRanges) {
        leftHeight += m_pimpl->stateVariables(static_cast<unsigned int>(range.offset) + 2);
    }
    leftHeight /= m_pimpl->leftPointRanges.size();

    double rightHeight = 0.0;
    for (auto& range : m_pimpl->rightPointRanges) {
        rightHeight += m_pimpl->stateVariables(static_cast<unsigned int>(range.offset) + 2);
    }
    rightHeight /= m_pimpl->rightPointRanges.size();

    constraint(0) = leftHeight - rightHeight;

    return true;
}

bool FeetRelativeHeightConstraint::constraintJacobianWRTState(double , const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->stateJacobianBuffer;

    return true;


}

bool FeetRelativeHeightConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t FeetRelativeHeightConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t FeetRelativeHeightConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool FeetRelativeHeightConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool FeetRelativeHeightConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}

bool FeetRelativeHeightConstraint::constraintSecondPartialDerivativeWRTState(double , const iDynTree::VectorDynSize &,
                                                                              const iDynTree::VectorDynSize &,
                                                                              const iDynTree::VectorDynSize &,
                                                                              iDynTree::MatrixDynSize &)
{
    return true;
}

bool FeetRelativeHeightConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                const iDynTree::VectorDynSize &/*control*/,
                                                                                const iDynTree::VectorDynSize &/*lambda*/,
                                                                                iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool FeetRelativeHeightConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                     const iDynTree::VectorDynSize &/*control*/,
                                                                                     const iDynTree::VectorDynSize &/*lambda*/,
                                                                                     iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool FeetRelativeHeightConstraint::constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool FeetRelativeHeightConstraint::constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool FeetRelativeHeightConstraint::constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}

