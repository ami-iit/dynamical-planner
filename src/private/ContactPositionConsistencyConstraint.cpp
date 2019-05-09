/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Constraints/ContactPositionConsistencyConstraint.h>
#include <DynamicalPlannerPrivate/Utilities/CheckEqualVector.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class ContactPositionConsistencyConstraint::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;

    std::string footName;
    iDynTree::FrameIndex footFrame;
    size_t contactIndex;
    iDynTree::Position positionInFoot;

    iDynTree::IndexRange positionPointRange;
    iDynTree::IndexRange jointsPositionRange, basePositionRange, baseQuaternionRange;

    iDynTree::VectorDynSize constraintValueBuffer;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::Position pointPosition;
    iDynTree::MatrixDynSize footJacobianBuffer, pointJacobianBuffer, stateJacobianBuffer, controlJacobianBuffer;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;
    iDynTree::MatrixFixSize<3, 6> /*footInternalTransformation,*/ footTransformationBuffer;

    iDynTree::Rotation baseRotation;
    iDynTree::Position basePosition;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    bool updateDoneOnceConstraint = false;
    bool updateDoneOnceStateJacobian = false;
    double tolerance;

    iDynTree::optimalcontrol::SparsityStructure stateJacobianSparsity, controlJacobianSparsity;
    iDynTree::optimalcontrol::SparsityStructure stateHessianSparsity, controlHessianSparsity, mixedHessianSparsity;

    levi::Expression asExpression, quaternionDerivative, jointsDerivative;
    std::vector<levi::Expression> quaternionQuaternionDerivatives, quaternionJointsDerivatives, jointsJointsDerivatives;
    iDynTree::VectorDynSize jointsHessianBuffer;

    void getRanges() {

        positionPointRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(contactIndex));
        assert(positionPointRange.isValid());

        basePositionRange = stateVariables.getIndexRange("BasePosition");
        assert(basePositionRange.isValid());

        baseQuaternionRange = stateVariables.getIndexRange("BaseQuaternion");
        assert(baseQuaternionRange.isValid());

        jointsPositionRange = stateVariables.getIndexRange("JointsPosition");
        assert(jointsPositionRange.isValid());
    }

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
        sharedKinDyn->updateRobotState(robotState);
    }

    void updateVariables (){
        updateRobotState();
        iDynTree::toEigen(pointPosition) = iDynTree::toEigen(stateVariables(positionPointRange));
    }

    bool sameVariables(bool updateDoneOnce) {
        bool same = updateDoneOnce;
        same = same && VectorsAreEqual(pointPosition, stateVariables(positionPointRange), tolerance);
        same = same && VectorsAreEqual(basePosition, stateVariables(basePositionRange), tolerance);
        same = same && VectorsAreEqual(baseQuaternion, stateVariables(baseQuaternionRange), tolerance);
        same = same && VectorsAreEqual(robotState.s, stateVariables(jointsPositionRange), tolerance);

        return same;
    }

    void setSparsity() {
        stateJacobianSparsity.clear();
        controlJacobianSparsity.clear();

        iDynTree::IndexRange fullRange;
        fullRange.offset = 0;
        fullRange.size = 3;

        stateJacobianSparsity.addIdentityBlock(0, static_cast<size_t>(basePositionRange.offset), 3);
        stateJacobianSparsity.addDenseBlock(fullRange, baseQuaternionRange);
        stateJacobianSparsity.addDenseBlock(fullRange, jointsPositionRange);
        stateJacobianSparsity.addIdentityBlock(0, static_cast<size_t>(positionPointRange.offset), 3);

        stateHessianSparsity.addDenseBlock(baseQuaternionRange, baseQuaternionRange);
        stateHessianSparsity.addDenseBlock(baseQuaternionRange, jointsPositionRange);
        stateHessianSparsity.addDenseBlock(jointsPositionRange, baseQuaternionRange);
        stateHessianSparsity.addDenseBlock(jointsPositionRange, jointsPositionRange);
        controlHessianSparsity.clear();
        mixedHessianSparsity.clear();

    }

    void clearDerivativesCache(std::vector<levi::Expression>& vector) {
        for (auto& expr : vector) {
            expr.clearDerivativesCache();
        }
    }

    ~Implementation() {
        asExpression.clearDerivativesCache();
        quaternionDerivative.clearDerivativesCache();
        jointsDerivative.clearDerivativesCache();
        clearDerivativesCache(quaternionQuaternionDerivatives);
        clearDerivativesCache(quaternionJointsDerivatives);
        clearDerivativesCache(jointsJointsDerivatives);
    }
};



ContactPositionConsistencyConstraint::ContactPositionConsistencyConstraint(const VariablesLabeller& stateVariables,
                                                                           const VariablesLabeller& controlVariables,
                                                                           std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                                                           std::shared_ptr<ExpressionsServer> expressionsServer,
                                                                           iDynTree::FrameIndex footFrame, const std::string &footName,
                                                                           const iDynTree::Position &positionInFoot, size_t contactIndex)
    : iDynTree::optimalcontrol::Constraint(3, "ContactPositionConsistency" + footName + std::to_string(contactIndex))
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());

    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;
    m_pimpl->footFrame = footFrame;
    assert(footFrame != iDynTree::FRAME_INVALID_INDEX);
    m_pimpl->footName = footName;
    m_pimpl->positionInFoot = positionInFoot;
    m_pimpl->contactIndex = contactIndex;

    m_pimpl->getRanges();

    m_pimpl->constraintValueBuffer.resize(3);
    m_pimpl->constraintValueBuffer.zero();
    m_pimpl->footJacobianBuffer.resize(6, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->pointJacobianBuffer.resize(3, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->stateJacobianBuffer.resize(3, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateJacobianBuffer.zero();
    m_pimpl->controlJacobianBuffer.resize(3, static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->tolerance = timelySharedKinDyn->getUpdateTolerance();

    iDynTree::toEigen(m_pimpl->footTransformationBuffer).leftCols<3>().setIdentity();


    m_isLowerBounded = true;
    m_isUpperBounded = true;
    m_upperBound.zero();
    m_lowerBound.zero();

    m_pimpl->setSparsity();

    m_pimpl->expressionsServer = expressionsServer;
    m_pimpl->jointsHessianBuffer.resize(static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));

    levi::Constant positionInFootExpr(iDynTree::toEigen(positionInFoot), footName + "_p_" + std::to_string(contactIndex));

    std::string frameName = timelySharedKinDyn->model().getFrameName(footFrame);
    m_pimpl->asExpression = m_pimpl->expressionsServer->worldToBase() *
        m_pimpl->expressionsServer->relativeTransform(timelySharedKinDyn->getFloatingBase(), frameName) * positionInFootExpr;

    m_pimpl->quaternionDerivative = m_pimpl->asExpression.getColumnDerivative(0, (m_pimpl->expressionsServer->baseQuaternion()));
    m_pimpl->jointsDerivative = m_pimpl->asExpression.getColumnDerivative(0, (m_pimpl->expressionsServer->jointsPosition()));

    for (Eigen::Index i = 0; i < 4; ++i) {

        m_pimpl->quaternionQuaternionDerivatives.push_back(m_pimpl->quaternionDerivative.getColumnDerivative(i, (m_pimpl->expressionsServer->baseQuaternion())));
        m_pimpl->quaternionJointsDerivatives.push_back(m_pimpl->quaternionDerivative.getColumnDerivative(i, (m_pimpl->expressionsServer->jointsPosition())));

    }

    for (Eigen::Index i = 0; i < m_pimpl->jointsPositionRange.size; ++i) {
        m_pimpl->jointsJointsDerivatives.push_back(m_pimpl->jointsDerivative.getColumnDerivative(i, (m_pimpl->expressionsServer->jointsPosition())));
    }

}

ContactPositionConsistencyConstraint::~ContactPositionConsistencyConstraint()
{ }

void ContactPositionConsistencyConstraint::setEqualityTolerance(double tolerance)
{
    assert(tolerance >= 0);

    iDynTree::toEigen(m_lowerBound).setConstant(-tolerance/2.0);
    iDynTree::toEigen(m_upperBound).setConstant(tolerance/2.0);
}

bool ContactPositionConsistencyConstraint::evaluateConstraint(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &constraint)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceConstraint))) {

//        m_pimpl->updateDoneOnceConstraint = true;
        m_pimpl->updateVariables();

        iDynTree::toEigen(m_pimpl->constraintValueBuffer) = iDynTree::toEigen(m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame) * m_pimpl->positionInFoot
                - m_pimpl->pointPosition);
    }

    constraint = m_pimpl->constraintValueBuffer;

    return true;
}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    if (!(m_pimpl->sameVariables(m_pimpl->updateDoneOnceStateJacobian))) {

//        m_pimpl->updateDoneOnceStateJacobian = true;
        m_pimpl->updateVariables();
        m_pimpl->expressionsServer->updateRobotState(time);

        bool ok = m_pimpl->sharedKinDyn->getFrameFreeFloatingJacobian(m_pimpl->robotState, m_pimpl->footFrame, m_pimpl->footJacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
        assert(ok);

        iDynTree::iDynTreeEigenMatrixMap footJacobianMap = iDynTree::toEigen(m_pimpl->footJacobianBuffer);
        iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);
        iDynTree::Transform footTransform = m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame);

//        Eigen::Map<Eigen::Matrix<double, 3, 6, Eigen::RowMajor> > footTransformationMap = iDynTree::toEigen(m_pimpl->footTransformationBuffer);

//        //    footTransformationMap = iDynTree::toEigen(m_pimpl->sharedKinDyn->getWorldTransform(m_pimpl->robotState, m_pimpl->footFrame).getRotation()) * iDynTree::toEigen(m_pimpl->footInternalTransformation);
//        footTransformationMap.rightCols<3>() = iDynTree::skew(iDynTree::toEigen(footTransform.getPosition() - (footTransform * m_pimpl->positionInFoot)));

//        iDynTree::iDynTreeEigenMatrixMap pointJacobianMap = iDynTree::toEigen(m_pimpl->pointJacobianBuffer);

//        //    iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));
//        iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap) = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

//        pointJacobianMap = footTransformationMap * footJacobianMap;

//        jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset) = pointJacobianMap.leftCols<3>();

//        jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) = pointJacobianMap.block<3, 3>(0, 3) * iDynTree::toEigen(m_pimpl->notNormalizedQuaternionMap);

//        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) = pointJacobianMap.rightCols(m_pimpl->jointsPositionRange.size);

//        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset).setIdentity();
//        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset) *= -1;


        iDynTree::Transform footTransformInBase = m_pimpl->sharedKinDyn->getBaseTransform(m_pimpl->robotState).inverse() * footTransform;
        iDynTree::Position expectedPointPositionInBase = footTransformInBase * m_pimpl->positionInFoot;
        iDynTree::Vector4 footQuaternion = footTransform.getRotation().asQuaternion();

        jacobianMap.block<3, 3>(0, m_pimpl->basePositionRange.offset).setIdentity();

        jacobianMap.block<3, 4>(0, m_pimpl->baseQuaternionRange.offset) =
                /*iDynTree::toEigen(RotatedVectorQuaternionJacobian(expectedPointPositionInBase, m_pimpl->baseQuaternionNormalized)) *
                iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion))*/ m_pimpl->quaternionDerivative.evaluate();

        jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size) =/*
                    footJacobianMap.topRightCorner(3, m_pimpl->jointsPositionRange.size) +
                    (iDynTree::toEigen(RotatedVectorQuaternionJacobian(m_pimpl->positionInFoot, footQuaternion)) *
                     iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivative(footQuaternion)) *
                     footJacobianMap.bottomRightCorner(3, m_pimpl->jointsPositionRange.size))*/ m_pimpl->jointsDerivative.evaluate();


        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset).setIdentity();
        jacobianMap.block<3,3>(0, m_pimpl->positionPointRange.offset) *= -1;
    }

    jacobian = m_pimpl->stateJacobianBuffer;

    return true;

}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::MatrixDynSize &jacobian)
{
    jacobian = m_pimpl->controlJacobianBuffer;
    return true;
}

size_t ContactPositionConsistencyConstraint::expectedStateSpaceSize() const
{
    return m_pimpl->stateVariables.size();
}

size_t ContactPositionConsistencyConstraint::expectedControlSpaceSize() const
{
    return m_pimpl->controlVariables.size();
}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateJacobianSparsity;
    return true;
}

bool ContactPositionConsistencyConstraint::constraintJacobianWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlJacobianSparsity;
    return true;
}

bool ContactPositionConsistencyConstraint::constraintSecondPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state,
                                                                                     const iDynTree::VectorDynSize &control,
                                                                                     const iDynTree::VectorDynSize &lambda,
                                                                                     iDynTree::MatrixDynSize &hessian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateVariables();
    m_pimpl->expressionsServer->updateRobotState(time);

    iDynTree::iDynTreeEigenMatrixMap hessianMap = iDynTree::toEigen(hessian);
    iDynTree::iDynTreeEigenConstVector lambdaMap = iDynTree::toEigen(lambda);

    iDynTree::iDynTreeEigenVector jointsMap = iDynTree::toEigen(m_pimpl->jointsHessianBuffer);
    Eigen::Matrix<double, 1, 4> quaternionHessian;

    for (Eigen::Index i = 0; i < 4; ++i) {

        quaternionHessian = lambdaMap.transpose() *
            m_pimpl->quaternionQuaternionDerivatives[static_cast<size_t>(i)].evaluate();

        hessianMap.block(m_pimpl->baseQuaternionRange.offset + i, m_pimpl->baseQuaternionRange.offset, 1, 4) = quaternionHessian;

        jointsMap =
            (m_pimpl->quaternionJointsDerivatives[static_cast<size_t>(i)].evaluate()).transpose() *
            lambdaMap;

        hessianMap.block(m_pimpl->baseQuaternionRange.offset + i, m_pimpl->jointsPositionRange.offset, 1, m_pimpl->jointsPositionRange.size) =
            jointsMap.transpose();

        hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->baseQuaternionRange.offset + i, m_pimpl->jointsPositionRange.size, 1) =
            jointsMap;
    }

    for (Eigen::Index i = 0; i < m_pimpl->jointsPositionRange.size; ++i) {
        jointsMap = (m_pimpl->jointsJointsDerivatives[static_cast<size_t>(i)].evaluate()).transpose() *
            lambdaMap;

        hessianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.offset + i, m_pimpl->jointsPositionRange.size, 1) =
            jointsMap;
    }

    return true;
}

bool ContactPositionConsistencyConstraint::constraintSecondPartialDerivativeWRTControl(double /*time*/, const iDynTree::VectorDynSize &/*state*/,
                                                                                       const iDynTree::VectorDynSize &/*control*/,
                                                                                       const iDynTree::VectorDynSize &/*lambda*/,
                                                                                       iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ContactPositionConsistencyConstraint::constraintSecondPartialDerivativeWRTStateControl(double /*time*/,
                                                                                            const iDynTree::VectorDynSize &/*state*/,
                                                                                            const iDynTree::VectorDynSize &/*control*/,
                                                                                            const iDynTree::VectorDynSize &/*lambda*/,
                                                                                            iDynTree::MatrixDynSize &/*hessian*/)
{
    return true;
}

bool ContactPositionConsistencyConstraint::constraintSecondPartialDerivativeWRTStateSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateHessianSparsity;
    return true;
}

bool ContactPositionConsistencyConstraint::constraintSecondPartialDerivativeWRTStateControlSparsity(iDynTree::optimalcontrol::SparsityStructure &stateControlSparsity)
{
    stateControlSparsity = m_pimpl->mixedHessianSparsity;
    return true;
}

bool ContactPositionConsistencyConstraint::constraintSecondPartialDerivativeWRTControlSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlHessianSparsity;
    return true;
}
