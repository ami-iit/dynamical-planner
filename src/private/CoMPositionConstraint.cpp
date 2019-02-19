/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Constraints/CoMPositionConstraint.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/CheckEqualVector.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/Core/MatrixDynSize.h>

#include <string>

#include <cassert>

using namespace DynamicalPlanner::Private;

class CoMPositionConstraint::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;
    iDynTree::IndexRange comPositionRange, jointsPositionRange, basePositionRange, baseQuaternionRange;
    iDynTree::VectorDynSize constraintValueBuffer;
    iDynTree::Position basePosition;
    iDynTree::Vector3 comPosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::Rotation baseRotation;
    iDynTree::MatrixDynSize comJacobianBuffer, stateJacobianBuffer, controlJacobianBuffer;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    bool updateDoneOnceConstraint = false;
    bool updateDoneOnceStateJacobian = false;
    double tolerance;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;
    levi::Expression asExpression, quaternionDerivative, jointsDerivative;
    iDynTree::VectorDynSize jointsHessianBuffer;

    void updateRobotState() {

        robotState = sharedKinDyn->currentState();

        iDynTree::toEigen(basePosition) = iDynTree::toEigen(stateVariables(basePositionRange));
        baseQuaternion = stateVariables(baseQuaternionRange);
        baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);

        robotState.base_quaternion = baseQuaternion;
        robotState.base_position = basePosition;

        robotState.s = stateVariables(jointsPositionRange);
    }

    void updateVariables (){
        updateRobotState();
        comPosition = stateVariables(comPositionRange);
    }

    bool sameVariables(bool updateDoneOnce) {
        bool same = updateDoneOnce;
        same = same && VectorsAreEqual(comPosition, stateVariables(comPositionRange), tolerance);
        same = same && VectorsAreEqual(basePosition, stateVariables(basePositionRange), tolerance);
        same = same && VectorsAreEqual(baseQuaternion, stateVariables(baseQuaternionRange), tolerance);
        same = same && VectorsAreEqual(robotState.s, stateVariables(jointsPositionRange), tolerance);
        return same;
    }

    void setSparsity() {
        stateSparsity.clear();
        controlSparsity.clear();

        iDynTree::IndexRange fullRange;
        fullRange.offset = 0;
        fullRange.size = 3;

        stateSparsity.addIdentityBlock(0, static_cast<size_t>(basePositionRange.offset), 3);
        stateSparsity.addDenseBlock(fullRange, baseQuaternionRange);
        stateSparsity.addDenseBlock(fullRange, jointsPositionRange);
        stateSparsity.addIdentityBlock(0, static_cast<size_t>(comPositionRange.offset), 3);
    }

};


CoMPositionConstraint::CoMPositionConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                             std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                             std::shared_ptr<ExpressionsServer> expressionsServer)
    : iDynTree::optimalcontrol::Constraint (3, "CoMPosition")
    , m_pimpl(std::make_unique<Implementation>())
{
    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());
    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;

    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    m_isLowerBounded = true;
    m_isUpperBounded = true;
    m_upperBound.zero();
    m_lowerBound.zero();

    m_pimpl->comPositionRange = m_pimpl->stateVariables.getIndexRange("CoMPosition");
    assert(m_pimpl->comPositionRange.isValid());

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->constraintValueBuffer.resize(3);
    m_pimpl->constraintValueBuffer.zero();
    m_pimpl->comJacobianBuffer.resize(3, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateJacobianBuffer.resize(3, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();
    m_pimpl->controlJacobianBuffer.resize(3, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->tolerance = timelySharedKinDyn->getUpdateTolerance();

    m_pimpl->setSparsity();

    m_pimpl->expressionsServer = expressionsServer;

    m_pimpl->asExpression = m_pimpl->expressionsServer->worldToBase() * m_pimpl->expressionsServer->comInBase();

    m_pimpl->quaternionDerivative = m_pimpl->asExpression.getColumnDerivative(0, m_pimpl->expressionsServer->baseQuaternion());
    m_pimpl->jointsDerivative = m_pimpl->asExpression.getColumnDerivative(0, m_pimpl->expressionsServer->jointsPosition());
    m_pimpl->jointsHessianBuffer.resize(static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
}

CoMPositionConstraint::~CoMPositionConstraint()
{ }

void CoMPositionConstraint::setEqualityTolerance(double tolerance)
{
    assert(tolerance >= 0);

    iDynTree::toEigen(m_lowerBound).setConstant(-tolerance/2.0);
    iDynTree::toEigen(m_upperBound).setConstant(tolerance/2.0);
}

bool CoMPositionConstraint::evaluateConstraint(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceConstraint))) {

//        m_pimpl->updateDoneOnceConstraint = true;
        m_pimpl->updateVariables();

        iDynTree::toEigen(m_pimpl->constraintValueBuffer) = iDynTree::toEigen(m_pimpl->sharedKinDyn->getCenterOfMassPosition(m_pimpl->robotState)) - iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->comPositionRange));

    }
    constraint = m_pimpl->constraintValueBuffer;

    return true;

}

bool CoMPositionConstraint::constraintJacobianWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceStateJacobian))) {

//        m_pimpl->updateDoneOnceStateJacobian = true;
        m_pimpl->updateVariables();

        bool ok = m_pimpl->sharedKinDyn->getCenterOfMassJacobian(m_pimpl->robotState, m_pimpl->comJacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
        assert(ok);

        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);
        iDynTree::iDynTreeEigenMatrixMap comJacobianMap = iDynTree::toEigen(m_pimpl->comJacobianBuffer);

        iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

        jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset) = comJacobianMap.topLeftCorner<3, 3>();
        jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) = comJacobianMap.block<3, 3>(0, 3) * iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap);
        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) = comJacobianMap.topRightCorner(3, m_pimpl->jointsPositionRange.size);

        jacobianMap.block<3,3>(0, m_pimpl->comPositionRange.offset).setIdentity();
        jacobianMap.block<3,3>(0, m_pimpl->comPositionRange.offset) *= -1;
    }

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;
}

bool CoMPositionConstraint::constraintJacobianWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/, const iDynTree::VectorDynSize &/*control*/, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t CoMPositionConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t CoMPositionConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool CoMPositionConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool CoMPositionConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}

bool CoMPositionConstraint::constraintSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                                      const iDynTree::VectorDynSize &control,
                                                                      const iDynTree::VectorDynSize &lambda, iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateVariables();
    m_pimpl->expressionsServer->updateRobotState(time, m_pimpl->robotState);

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(hessian);
    iDynTree::iDynTreeEigenConstVector lambdaMap = iDynTree::toEigen(lambda);

    iDynTree::iDynTreeEigenVector jointsMap = iDynTree::toEigen(m_pimpl->jointsHessianBuffer);
    Eigen::Matrix<double, 1, 4> quaternionHessian;

    for (Eigen::Index i = 0; i < 4; ++i) {

        quaternionHessian = lambdaMap.transpose() *
            m_pimpl->quaternionDerivative.getColumnDerivative(i, (m_pimpl->expressionsServer->baseQuaternion())).evaluate();

        hessianMap.block(m_pimpl->baseQuaternionRange.offset + i, m_pimpl->baseQuaternionRange.offset, 1, 4) = quaternionHessian;

        jointsMap =
            (m_pimpl->quaternionDerivative.getColumnDerivative(i, (m_pimpl->expressionsServer->jointsPosition())).evaluate()).transpose() *
            lambdaMap;

        hessianMap.block(m_pimpl->baseQuaternionRange.offset + i, m_pimpl->jointsPositionRange.offset, 1, m_pimpl->jointsPositionRange.size) =
            jointsMap.transpose();

        hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionRange.offset + i, m_pimpl->jointsPositionRange.size, 1) =
            jointsMap;
    }

    for (Eigen::Index i = 0; i < m_pimpl->jointsPositionRange.size; ++i) {
        jointsMap = (m_pimpl->jointsDerivative.getColumnDerivative(i, (m_pimpl->expressionsServer->jointsPosition())).evaluate()).transpose() *
            lambdaMap;

        hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.offset + i, m_pimpl->jointsPositionRange.size, 1) =
            jointsMap;
    }

    return true;
}

bool CoMPositionConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                        const iDynTree::VectorDynSize &/*control*/,
                                                                        const iDynTree::VectorDynSize &/*lambda*/,
                                                                        iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool CoMPositionConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                             const iDynTree::VectorDynSize &/*control*/,
                                                                             const iDynTree::VectorDynSize &/*lambda*/,
                                                                             iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}
