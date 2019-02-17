/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Constraints/FeetLateralDistanceConstraint.h>
#include <DynamicalPlannerPrivate/Utilities/CheckEqualVector.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class FeetLateralDistanceConstraint::Implementation {
public:
    VariablesLabeller stateVariables, controlVariables;

    unsigned int lateralIndex;
    iDynTree::FrameIndex referenceFootFrame, otherFootFrame;
    iDynTree::IndexRange jointsPositionRange;

    iDynTree::MatrixDynSize relativeJacobianBuffer, stateJacobianBuffer, controlJacobianBuffer;
    double feetDistance = 0.0;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    bool updateDoneOnceConstraint = false;
    bool updateDoneOnceStateJacobian = false;
    double tolerance;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;

    levi::Expression asExpression, jointsDerivative;
    iDynTree::VectorDynSize jointsHessianBuffer;

    void updateRobotState() {

        robotState = sharedKinDyn->currentState();

        robotState.s = stateVariables(jointsPositionRange);
    }

    bool sameVariables(bool updateDoneOnce) {
        bool same = updateDoneOnce;
        same = same && VectorsAreEqual(robotState.s, stateVariables(jointsPositionRange), tolerance);
        return same;
    }

    void setSparsity() {
        stateSparsity.clear();
        controlSparsity.clear();

        iDynTree::IndexRange fullRange;
        fullRange.offset = 0;
        fullRange.size = 1;

        stateSparsity.addDenseBlock(fullRange, jointsPositionRange);
    }

};


FeetLateralDistanceConstraint::FeetLateralDistanceConstraint(const VariablesLabeller& stateVariables,
                                                             const VariablesLabeller& controlVariables,
                                                             std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                                             std::shared_ptr<ExpressionsServer> expressionsServer,
                                                             unsigned int lateralIndex, iDynTree::FrameIndex referenceFootFrame,
                                                             iDynTree::FrameIndex otherFootFrame)
    : iDynTree::optimalcontrol::Constraint (1, "FeetLateralConstraint")
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());
    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;
    assert(lateralIndex < 3);
    m_pimpl->lateralIndex = lateralIndex;
    assert(referenceFootFrame != iDynTree::FRAME_INVALID_INDEX);
    assert(otherFootFrame != iDynTree::FRAME_INVALID_INDEX);
    m_pimpl->referenceFootFrame = referenceFootFrame;
    m_pimpl->otherFootFrame = otherFootFrame;

    m_pimpl->jointsPositionRange = stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->relativeJacobianBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateJacobianBuffer.resize(1, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();

    m_pimpl->controlJacobianBuffer.resize(1, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->tolerance = timelySharedKinDyn->getUpdateTolerance();

    m_lowerBound(0) = 0.1;

    m_isLowerBounded = true;
    m_isUpperBounded = false;
    m_lowerBound.zero();

    m_pimpl->setSparsity();

    m_pimpl->expressionsServer = expressionsServer;
    m_pimpl->jointsHessianBuffer.resize(static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));

    std::string referenceFrameName = timelySharedKinDyn->model().getFrameName(referenceFootFrame);
    std::string otherFootName = timelySharedKinDyn->model().getFrameName(otherFootFrame);

    m_pimpl->asExpression = (m_pimpl->expressionsServer->relativePosition(referenceFrameName, otherFootName))->row(lateralIndex);
    m_pimpl->jointsDerivative = m_pimpl->asExpression.getColumnDerivative(0, *(m_pimpl->expressionsServer->jointsPosition()));

}

FeetLateralDistanceConstraint::~FeetLateralDistanceConstraint()
{ }

bool FeetLateralDistanceConstraint::setMinimumDistance(double minDistance)
{
    if (minDistance < 0)
        return false;

    m_lowerBound(0) = minDistance;

    return true;
}

bool FeetLateralDistanceConstraint::evaluateConstraint(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);
    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceConstraint))) {

//        m_pimpl->updateDoneOnceConstraint = true;
        m_pimpl->updateRobotState();

        m_pimpl->feetDistance = m_pimpl->sharedKinDyn->getRelativeTransform(m_pimpl->robotState,
                                                                   m_pimpl->referenceFootFrame,
                                                                   m_pimpl->otherFootFrame).getPosition()(m_pimpl->lateralIndex);
    }
    constraint(0) = m_pimpl->feetDistance;
    return true;
}

bool FeetLateralDistanceConstraint::constraintJacobianWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceStateJacobian))) {

//        m_pimpl->updateDoneOnceStateJacobian = true;
        m_pimpl->updateRobotState();

        bool ok = m_pimpl->sharedKinDyn->getRelativeJacobian(m_pimpl->robotState, m_pimpl->referenceFootFrame,
                                                             m_pimpl->otherFootFrame, m_pimpl->relativeJacobianBuffer,
                                                             iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
        assert(ok);

        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);
        iDynTree::iDynTreeEigenMatrixMap relativeJacobianMap = iDynTree::toEigen(m_pimpl->relativeJacobianBuffer);

        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 1, m_pimpl->jointsPositionRange.size) =
                relativeJacobianMap.block(m_pimpl->lateralIndex, 0, 1, m_pimpl->jointsPositionRange.size);
    }

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;


}

bool FeetLateralDistanceConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t FeetLateralDistanceConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t FeetLateralDistanceConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool FeetLateralDistanceConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool FeetLateralDistanceConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}

bool FeetLateralDistanceConstraint::constraintSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                                              const iDynTree::VectorDynSize &control,
                                                                              const iDynTree::VectorDynSize &lambda,
                                                                              iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
    m_pimpl->expressionsServer->updateRobotState(time, m_pimpl->robotState);

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(hessian);
    iDynTree::iDynTreeEigenConstVector lambdaMap = iDynTree::toEigen(lambda);

    iDynTree::iDynTreeEigenVector jointsMap = iDynTree::toEigen(m_pimpl->jointsHessianBuffer);
    Eigen::Matrix<double, 1, 4> quaternionHessian;

    for (Eigen::Index i = 0; i < m_pimpl->jointsPositionRange.size; ++i) {
        jointsMap = (m_pimpl->jointsDerivative.getColumnDerivative(i, *(m_pimpl->expressionsServer->jointsPosition())).evaluate()).transpose() *
            lambdaMap;

        hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.offset + i, m_pimpl->jointsPositionRange.size, 1) =
            jointsMap;
    }

    return true;
}

bool FeetLateralDistanceConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                const iDynTree::VectorDynSize &/*control*/,
                                                                                const iDynTree::VectorDynSize &/*lambda*/,
                                                                                iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool FeetLateralDistanceConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                     const iDynTree::VectorDynSize &/*control*/,
                                                                                     const iDynTree::VectorDynSize &/*lambda*/,
                                                                                     iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}
