/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Constraints/CentroidalMomentumConstraint.h>
#include <DynamicalPlannerPrivate/Utilities/CheckEqualVector.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <DynamicalPlannerPrivate/Utilities/levi/CoMInBaseExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeVelocityExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/MomentumInBaseExpression.h>
#include <cassert>
#include <thread>
#include <future>
#include <chrono>

using namespace DynamicalPlanner::Private;

class CentroidalMomentumConstraint::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;

    iDynTree::IndexRange momentumRange, comPositionRange, basePositionRange, baseQuaternionRange, jointsPositionRange, jointsVelocityRange;
//    iDynTree::IndexRange baseVelocityRange;
    iDynTree::IndexRange baseLinearVelocityRange, baseQuaternionDerivativeRange;
    iDynTree::VectorDynSize constraintValueBuffer;
    iDynTree::Position basePosition;
    iDynTree::Position comPosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized, baseQuaternionVelocity;
    iDynTree::Rotation baseRotation;
    iDynTree::Transform comTransform;
    iDynTree::Vector6 momentum;

    iDynTree::MatrixDynSize cmmMatrixInCoMBuffer, cmmMatrixInBaseBuffer,
    momentumDerivativeBuffer, stateJacobianBuffer, controlJacobianBuffer, comJacobianBuffer;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    bool updateDoneOnceConstraint = false;
    bool updateDoneOnceStateJacobian = false;
    bool updateDoneOnceControlJacobian = false;
    double tolerance;

    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;

    levi::Expression asExpression, quaternionDerivative, jointsDerivative;

    levi::Variable normalizedQuaternion;

    levi::Variable comPositionVariable;

    levi::Variable lagrangeMultipliers = levi::Variable(3, "lambdaCentroidal");
    levi::Expression quatQuatHessian, quatJointsHessian, jointsJointsHessian, quatLinVelHessian, quatQuatVelHessian, quatJointsVelHessian,
        jointsLinVelHessian, jointsQuatVelHessian, jointsJointsVelHessian;


    void getRanges() {

        momentumRange = stateVariables.getIndexRange("Momentum");
        assert(momentumRange.isValid());

        comPositionRange = stateVariables.getIndexRange("CoMPosition");
        assert(comPositionRange.isValid());

        basePositionRange = stateVariables.getIndexRange("BasePosition");
        assert(basePositionRange.isValid());

        baseQuaternionRange = stateVariables.getIndexRange("BaseQuaternion");
        assert(baseQuaternionRange.isValid());

        jointsPositionRange = stateVariables.getIndexRange("JointsPosition");
        assert(jointsPositionRange.isValid());

//        baseVelocityRange = controlVariables.getIndexRange("BaseVelocity");
//        assert(baseVelocityRange.isValid());

        baseLinearVelocityRange = controlVariables.getIndexRange("BaseLinearVelocity");
        assert(baseLinearVelocityRange.isValid());

        baseQuaternionDerivativeRange = controlVariables.getIndexRange("BaseQuaternionDerivative");
        assert(baseQuaternionDerivativeRange.isValid());

        jointsVelocityRange = controlVariables.getIndexRange("JointsVelocity");
        assert(jointsVelocityRange.isValid());
    }

    void updateRobotState() {

        robotState = sharedKinDyn->currentState();

        iDynTree::toEigen(basePosition) = iDynTree::toEigen(stateVariables(basePositionRange));
        baseQuaternion = stateVariables(baseQuaternionRange);
        baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);
        iDynTree::toEigen(baseQuaternionVelocity) = iDynTree::toEigen(controlVariables(baseQuaternionDerivativeRange));

        robotState.base_quaternion = baseQuaternion;
        robotState.base_position = basePosition;

        robotState.s = stateVariables(jointsPositionRange);

        robotState.s_dot = controlVariables(jointsVelocityRange);


        robotState.base_linearVelocity = controlVariables(baseLinearVelocityRange);
        robotState.base_quaternionVelocity = baseQuaternionVelocity;

    }

    void updateVariables (){
        updateRobotState();
//        comPosition = stateVariables(comPositionRange);
//        iDynTree::toEigen(comPositionInverse) = -1 * iDynTree::toEigen(comPosition);
        comPosition = sharedKinDyn->getCenterOfMassPosition(robotState);
        comTransform.setPosition(iDynTree::Position(0.0, 0.0, 0.0) - comPosition);
        comTransform.setRotation(iDynTree::Rotation::Identity());
        momentum = stateVariables(momentumRange);
    }

    bool sameVariables(bool updateDoneOnce) {
        bool same = updateDoneOnce;
        same = same && VectorsAreEqual(momentum, stateVariables(momentumRange), tolerance);
        same = same && VectorsAreEqual(comPosition, stateVariables(comPositionRange), tolerance);
        same = same && VectorsAreEqual(basePosition, stateVariables(basePositionRange), tolerance);
        same = same && VectorsAreEqual(baseQuaternion, stateVariables(baseQuaternionRange), tolerance);
        same = same && VectorsAreEqual(robotState.s, stateVariables(jointsPositionRange), tolerance);
        same = same && VectorsAreEqual(robotState.base_linearVelocity, controlVariables(baseLinearVelocityRange), tolerance);
        same = same && VectorsAreEqual(baseQuaternionVelocity, controlVariables(baseQuaternionDerivativeRange), tolerance);
        same = same && VectorsAreEqual(robotState.s_dot, controlVariables(jointsVelocityRange), tolerance);

        return same;
    }

    void setSparsity() {
        stateSparsity.clear();
        controlSparsity.clear();

        iDynTree::IndexRange fullRange;
        fullRange.size = 3;
        fullRange.offset = 0;

        stateSparsity.addDenseBlock(fullRange, jointsPositionRange);
        stateSparsity.addIdentityBlock(0, static_cast<size_t>(momentumRange.offset) + 3, 3);
        stateSparsity.addDenseBlock(fullRange, baseQuaternionRange);

        controlSparsity.addDenseBlock(fullRange, baseQuaternionDerivativeRange);
        controlSparsity.addDenseBlock(fullRange, baseLinearVelocityRange);
        controlSparsity.addDenseBlock(fullRange, jointsVelocityRange);

    }

    void constructExpressions() {
        const iDynTree::Model& model = timedSharedKinDyn->model();
        std::string baseFrame = timedSharedKinDyn->getFloatingBase();
        iDynTree::LinkIndex baseIndex = model.getLinkIndex(baseFrame);
        assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

        iDynTree::LinkConstPtr baseLink = model.getLink(baseIndex);
        levi::Constant baseInertia(iDynTree::toEigen(baseLink->getInertia().asMatrix()),"I_b");

        normalizedQuaternion = expressionsServer->normalizedBaseQuaternion().asIndependentVariable();
        comPositionVariable = expressionsServer->comInBase().asIndependentVariable();

        levi::Expression baseTwist = BodyTwistFromQuaternionVelocity(expressionsServer->baseLinearVelocity(),
                                                                     expressionsServer->baseQuaternionVelocity(),
                                                                     normalizedQuaternion, "baseTwist");

        levi::Expression worldToBaseRotation = RotationExpression(normalizedQuaternion);

        levi::Expression comInBasePosition = comPositionVariable;

        levi::Expression mixedAdjointBottomRows = worldToBaseRotation *
            levi::Expression::Horzcat((-comInBasePosition).skew(), levi::Identity(3,3), "G[b]_X_b");

        asExpression = mixedAdjointBottomRows * DynamicalPlanner::Private::MomentumInBaseExpression(expressionsServer.get(), baseTwist.asVariable());

        levi::Expression notNormalizedQuaternionMapExpr =
            expressionsServer->normalizedBaseQuaternion().getColumnDerivative(0, expressionsServer->baseQuaternion());

        quaternionDerivative = asExpression.getColumnDerivative(0, normalizedQuaternion) * notNormalizedQuaternionMapExpr;

        levi::Expression comJacobian = expressionsServer->comInBase().getColumnDerivative(0, expressionsServer->jointsPosition());

        jointsDerivative = asExpression.getColumnDerivative(0, expressionsServer->jointsPosition())
            + asExpression.getColumnDerivative(0, comPositionVariable) * comJacobian;

        levi::Expression lagrangian = lagrangeMultipliers.transpose() * asExpression;
        levi::Expression jointsLagrangian = (lagrangian.getColumnDerivative(0, expressionsServer->jointsPosition())
                                             + lagrangian.getColumnDerivative(0, comPositionVariable) * comJacobian).transpose();
        levi::Expression quaternionLagrangian = (lagrangian.getColumnDerivative(0, normalizedQuaternion) * notNormalizedQuaternionMapExpr).transpose();

        quatQuatHessian = quaternionLagrangian.getColumnDerivative(0, expressionsServer->baseQuaternion()) +
            quaternionLagrangian.getColumnDerivative(0, normalizedQuaternion) * notNormalizedQuaternionMapExpr;
        quatJointsHessian = quaternionLagrangian.getColumnDerivative(0, expressionsServer->jointsPosition()) +
            quaternionLagrangian.getColumnDerivative(0, comPositionVariable) * comJacobian;

        jointsJointsHessian = jointsLagrangian.getColumnDerivative(0, expressionsServer->jointsPosition()) +
            jointsLagrangian.getColumnDerivative(0, comPositionVariable) * comJacobian;

        quatLinVelHessian = quaternionLagrangian.getColumnDerivative(0, expressionsServer->baseLinearVelocity());
        quatQuatVelHessian = quaternionLagrangian.getColumnDerivative(0, expressionsServer->baseQuaternionVelocity());
        quatJointsVelHessian = quaternionLagrangian.getColumnDerivative(0, expressionsServer->jointsVelocity());

        jointsLinVelHessian = jointsLagrangian.getColumnDerivative(0, expressionsServer->baseLinearVelocity());
        jointsQuatVelHessian = jointsLagrangian.getColumnDerivative(0, expressionsServer->baseQuaternionVelocity());
        jointsJointsVelHessian = jointsLagrangian.getColumnDerivative(0, expressionsServer->jointsVelocity());
    }

    ~Implementation() {
        asExpression.clearDerivativesCache();
        quaternionDerivative.clearDerivativesCache();
        jointsDerivative.clearDerivativesCache();
        quatQuatHessian.clearDerivativesCache();
        quatJointsHessian.clearDerivativesCache();
        jointsJointsHessian.clearDerivativesCache();
        quatLinVelHessian.clearDerivativesCache();
        quatQuatVelHessian.clearDerivativesCache();
        quatJointsVelHessian.clearDerivativesCache();
        jointsLinVelHessian.clearDerivativesCache();
        jointsQuatVelHessian.clearDerivativesCache();
        jointsJointsVelHessian.clearDerivativesCache();
    }

};



CentroidalMomentumConstraint::CentroidalMomentumConstraint(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables,
                                                           std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                                           std::shared_ptr<ExpressionsServer> expressionServer)
    : iDynTree::optimalcontrol::Constraint (3, "CentroidalMomentum")
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

    m_pimpl->getRanges();

    m_pimpl->constraintValueBuffer.resize(6);
    m_pimpl->constraintValueBuffer.zero();
    m_pimpl->cmmMatrixInCoMBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->cmmMatrixInCoMBuffer.zero();
    m_pimpl->cmmMatrixInBaseBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->cmmMatrixInBaseBuffer.zero();
    m_pimpl->momentumDerivativeBuffer.resize(6, static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateJacobianBuffer.resize(6, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();
    m_pimpl->controlJacobianBuffer.resize(6, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();
    m_pimpl->comJacobianBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->comJacobianBuffer.zero();

    m_pimpl->tolerance = timelySharedKinDyn->getUpdateTolerance();

    m_pimpl->setSparsity();

    m_pimpl->expressionsServer = expressionServer;

    m_pimpl->constructExpressions();

}

void CentroidalMomentumConstraint::setEqualityTolerance(double tolerance)
{
    assert(tolerance >= 0);

    iDynTree::toEigen(m_lowerBound).setConstant(-tolerance/2.0);
    iDynTree::toEigen(m_upperBound).setConstant(tolerance/2.0);

}

CentroidalMomentumConstraint::~CentroidalMomentumConstraint()
{ }

bool CentroidalMomentumConstraint::evaluateConstraint(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceConstraint))) {

//        m_pimpl->updateDoneOnceConstraint = true;
        m_pimpl->updateVariables();

        iDynTree::SpatialMomentum expectedMomentum;
        expectedMomentum = m_pimpl->comTransform *
                (m_pimpl->sharedKinDyn->getBaseTransform(m_pimpl->robotState) *
                 m_pimpl->sharedKinDyn->getLinearAngularMomentum(m_pimpl->robotState, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION));

        iDynTree::toEigen(m_pimpl->constraintValueBuffer) = iDynTree::toEigen(expectedMomentum) - iDynTree::toEigen(m_pimpl->momentum);
    }

    iDynTree::toEigen(constraint) = iDynTree::toEigen(m_pimpl->constraintValueBuffer).bottomRows<3>();

    return true;
}

bool CentroidalMomentumConstraint::constraintJacobianWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceStateJacobian))) {

//        m_pimpl->updateDoneOnceStateJacobian = true;
        m_pimpl->updateVariables();

        iDynTree::Transform G_T_B = m_pimpl->comTransform * m_pimpl->sharedKinDyn->getBaseTransform(m_pimpl->robotState);

        iDynTree::SpatialMomentum momentumInCoM, momentumInBase;
        momentumInBase = m_pimpl->sharedKinDyn->getLinearAngularMomentum(m_pimpl->robotState, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        momentumInCoM = G_T_B * momentumInBase;


        bool ok = m_pimpl->sharedKinDyn->getLinearAngularMomentumJointsDerivative(m_pimpl->robotState, m_pimpl->momentumDerivativeBuffer);
        assert(ok);

        ok = m_pimpl->sharedKinDyn->getCenterOfMassJacobian(m_pimpl->robotState, m_pimpl->comJacobianBuffer,
                                                                 iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
        assert(ok);

        iDynTree::iDynTreeEigenMatrixMap comJacobianMap = iDynTree::toEigen(m_pimpl->comJacobianBuffer);

        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);

        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> skewMomentum = iDynTree::skew(iDynTree::toEigen(momentumInCoM).topRows<3>());

        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 6, m_pimpl->jointsPositionRange.size) =
                iDynTree::toEigen(G_T_B.asAdjointTransformWrench()) * iDynTree::toEigen(m_pimpl->momentumDerivativeBuffer);


        jacobianMap.block(3, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) +=
                skewMomentum * comJacobianMap.rightCols(m_pimpl->jointsPositionRange.size);

        jacobianMap.block<6,6>(0, m_pimpl->momentumRange.offset).setIdentity();

        jacobianMap.block<6,6>(0, m_pimpl->momentumRange.offset) *= -1;

//        jacobianMap.block<3,3>(3, m_pimpl->comPositionRange.offset) = iDynTree::skew(iDynTree::toEigen(momentumInCoM).topRows<3>());

        iDynTree::Matrix4x4 normalizedQuaternionDerivative = NormalizedQuaternionDerivative(m_pimpl->baseQuaternion);

        iDynTree::MatrixFixSize<3,4> linearPartDerivative;

        iDynTree::toEigen(linearPartDerivative) = iDynTree::toEigen(RotatedVectorQuaternionJacobian(momentumInBase.getLinearVec3(), m_pimpl->baseQuaternionNormalized)) *
                iDynTree::toEigen(normalizedQuaternionDerivative);

        jacobianMap.block<3,4>(0, m_pimpl->baseQuaternionRange.offset) = iDynTree::toEigen(linearPartDerivative);

        iDynTree::Position comPositionInBase = m_pimpl->sharedKinDyn->getBaseTransform(m_pimpl->robotState).inverse() * m_pimpl->comPosition;
        iDynTree::Vector3 comCrossMomentum;

        iDynTree::toEigen(comCrossMomentum) = (- iDynTree::toEigen(comPositionInBase)).cross(iDynTree::toEigen(momentumInBase.getLinearVec3()));

        jacobianMap.block<3,4>(3, m_pimpl->baseQuaternionRange.offset) = (iDynTree::toEigen(RotatedVectorQuaternionJacobian(momentumInBase.getAngularVec3(), m_pimpl->baseQuaternionNormalized))
                                                                           + iDynTree::toEigen(RotatedVectorQuaternionJacobian(comCrossMomentum, m_pimpl->baseQuaternionNormalized))) * iDynTree::toEigen(normalizedQuaternionDerivative);


         ok = m_pimpl->sharedKinDyn->getLinearAngularMomentumJacobian(m_pimpl->robotState, m_pimpl->cmmMatrixInBaseBuffer, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        assert(ok);

        jacobianMap.block<6,4>(0, m_pimpl->baseQuaternionRange.offset) += iDynTree::toEigen(G_T_B.asAdjointTransformWrench()) *
                (iDynTree::toEigen(m_pimpl->cmmMatrixInBaseBuffer).block<6,3>(0,3) * iDynTree::toEigen(
                     QuaternionLeftTrivializedDerivativeInverseTimesQuaternionDerivativeJacobian(m_pimpl->baseQuaternionVelocity)) *
                 iDynTree::toEigen(normalizedQuaternionDerivative));


//        m_pimpl->expressionsServer->updateRobotState(time);
//        jacobianMap.block<3,4>(3, m_pimpl->baseQuaternionRange.offset) = m_pimpl->quaternionDerivative.evaluate();
//        jacobianMap.block(3, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) = m_pimpl->jointsDerivative.evaluate();

    }

    iDynTree::toEigen(jacobian) = iDynTree::toEigen(m_pimpl->stateJacobianBuffer).bottomRows<3>();

    return true;
}

bool CentroidalMomentumConstraint::constraintJacobianWRTControl(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceControlJacobian))) {

//        m_pimpl->updateDoneOnceControlJacobian = true;
        m_pimpl->updateVariables();

        iDynTree::Transform G_T_B = m_pimpl->comTransform * m_pimpl->sharedKinDyn->getBaseTransform(m_pimpl->robotState);

        bool ok = m_pimpl->sharedKinDyn->getLinearAngularMomentumJacobian(m_pimpl->robotState, m_pimpl->cmmMatrixInBaseBuffer, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        assert(ok);

        iDynTree::toEigen(m_pimpl->cmmMatrixInCoMBuffer) = iDynTree::toEigen(G_T_B.asAdjointTransformWrench()) * iDynTree::toEigen(m_pimpl->cmmMatrixInBaseBuffer);

        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->controlJacobianBuffer);

        jacobianMap.block<6,3>(0, m_pimpl->baseLinearVelocityRange.offset) = iDynTree::toEigen(m_pimpl->cmmMatrixInCoMBuffer).leftCols<3>();

        jacobianMap.block<6,4>(0, m_pimpl->baseQuaternionDerivativeRange.offset) = iDynTree::toEigen(m_pimpl->cmmMatrixInCoMBuffer).block<6,3>(0, 3) *
                iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized));

        jacobianMap.block(0, m_pimpl->jointsVelocityRange.offset, 6, m_pimpl->jointsVelocityRange.size) = iDynTree::toEigen(m_pimpl->cmmMatrixInCoMBuffer).rightCols(m_pimpl->jointsVelocityRange.size);
    }
    iDynTree::toEigen(jacobian) = iDynTree::toEigen(m_pimpl->controlJacobianBuffer).bottomRows<3>();

    return true;
}

size_t CentroidalMomentumConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t CentroidalMomentumConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool CentroidalMomentumConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool CentroidalMomentumConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}

bool CentroidalMomentumConstraint::constraintSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, const iDynTree::VectorDynSize &lambda, iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateVariables();
    m_pimpl->expressionsServer->updateRobotState(time);

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(hessian);
    iDynTree::iDynTreeEigenConstVector lambdaMap = iDynTree::toEigen(lambda);

    m_pimpl->lagrangeMultipliers = lambdaMap;

    auto quatQuat = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->quatQuatHessian);
    auto quatJoints = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->quatJointsHessian);
    auto jointsJoints = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->jointsJointsHessian);

    using namespace std::chrono_literals;
    std::this_thread::sleep_for(10us);

    bool quatQuatDone = false, quatJointsDone = false, jointsJointsDone = false;

    while (!quatQuatDone || !quatJointsDone || !jointsJointsDone) {

        std::this_thread::sleep_for(1us);

        if (!quatQuatDone) {
            if (quatQuat.wait_for(1us) == std::future_status::ready) {
                hessianMap.block<4,4>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseQuaternionRange.offset) = quatQuat.get();
                quatQuatDone = true;
            }
        }

        if (!quatJointsDone) {
            if (quatJoints.wait_for(1us) == std::future_status::ready) {
                const Eigen::MatrixXd hessianOutput = quatJoints.get();
                hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionRange.offset, m_pimpl->jointsPositionRange.size, 4) =
                    hessianOutput.transpose();
                hessianMap.block(m_pimpl->baseQuaternionRange.offset, m_pimpl->jointsPositionRange.offset, 4, m_pimpl->jointsPositionRange.size) =
                    hessianOutput;
                quatJointsDone = true;
            }
        }

        if (!jointsJointsDone) {
            if (jointsJoints.wait_for(1us) == std::future_status::ready) {
                hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.offset,
                                 m_pimpl->jointsPositionRange.size, m_pimpl->jointsPositionRange.size) = jointsJoints.get();
                jointsJointsDone = true;
            }
        }
    }

    return true;
}

bool CentroidalMomentumConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                               const iDynTree::VectorDynSize &/*control*/,
                                                                               const iDynTree::VectorDynSize &/*lambda*/,
                                                                               iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool CentroidalMomentumConstraint::constraintSecondPartialDerivativeWRTStateControl(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, const iDynTree::VectorDynSize &lambda, iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateVariables();
    m_pimpl->expressionsServer->updateRobotState(time);

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(hessian);
    iDynTree::iDynTreeEigenConstVector lambdaMap = iDynTree::toEigen(lambda);

    m_pimpl->lagrangeMultipliers = lambdaMap;

    auto quatLinVel = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->quatLinVelHessian);
    auto quatQuatVel = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->quatQuatVelHessian);
    auto quatJointsVel = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->quatJointsVelHessian);
    auto jointsLinVel = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->jointsLinVelHessian);
    auto jointsQuatVel = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->jointsQuatVelHessian);
    auto jointsJointsVel = std::async(std::launch::async, &levi::Expression::evaluate, m_pimpl->jointsJointsVelHessian);

    using namespace std::chrono_literals;
    std::this_thread::sleep_for(10us);

    bool quatLinVelDone = false;
    bool quatQuatVelDone = false;
    bool quatJointsVelDone = false;
    bool jointsLinVelDone = false;
    bool jointsQuatVelDone = false;
    bool jointsJointsVelDone = false;

    while (!quatLinVelDone ||
           !quatQuatVelDone ||
           !quatJointsVelDone ||
           !jointsLinVelDone ||
           !jointsQuatVelDone ||
           !jointsJointsVelDone) {

        std::this_thread::sleep_for(1us);

        if (!quatLinVelDone) {
            if (quatLinVel.wait_for(1us) == std::future_status::ready) {
                hessianMap.block<4, 3>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseLinearVelocityRange.offset) = quatLinVel.get();
                quatLinVelDone = true;
            }
        }

        if (!quatQuatVelDone) {
            if (quatQuatVel.wait_for(1us) == std::future_status::ready) {
                hessianMap.block<4,4>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseQuaternionDerivativeRange.offset) = quatQuatVel.get();
                quatQuatVelDone = true;
            }
        }

        if (!quatJointsVelDone) {
            if (quatJointsVel.wait_for(1us) == std::future_status::ready) {
                hessianMap.block(m_pimpl->baseQuaternionRange.offset, m_pimpl->jointsVelocityRange.offset,
                                 4, m_pimpl->jointsVelocityRange.size) = quatJointsVel.get();
                quatJointsVelDone = true;
            }
        }

        if (!jointsLinVelDone) {
            if (jointsLinVel.wait_for(1us) == std::future_status::ready) {
                hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseLinearVelocityRange.offset,
                                 m_pimpl->jointsPositionRange.size, 3) = jointsLinVel.get();
                jointsLinVelDone = true;
            }
        }

        if (!jointsQuatVelDone) {
            if (jointsQuatVel.wait_for(1us) == std::future_status::ready) {
                hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionDerivativeRange.offset,
                                 m_pimpl->jointsPositionRange.size, 4) = jointsQuatVel.get();
                jointsQuatVelDone = true;
            }
        }

        if (!jointsJointsVelDone) {
            if (jointsJointsVel.wait_for(1us) == std::future_status::ready) {
                hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsVelocityRange.offset,
                                 m_pimpl->jointsPositionRange.size, m_pimpl->jointsVelocityRange.size) = jointsJointsVel.get();
                jointsJointsVelDone = true;
            }
        }
    }

    return true;
}
