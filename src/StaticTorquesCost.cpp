/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/StaticTorquesCost.h>
#include <DynamicalPlannerPrivate/QuaternionUtils.h>

#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class StaticTorquesCost::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;

    iDynTree::Position basePosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized;
    iDynTree::Rotation baseRotation;
    iDynTree::MatrixFixSize<3, 4> notNormalizedQuaternionMap;

    iDynTree::VectorDynSize stateGradientBuffer, controlGradientBuffer;
    iDynTree::MatrixDynSize fullJacobianBuffer, jointsJacobianBuffer, massMatrixBuffer;

    typedef struct {
        iDynTree::FrameIndex footFrame;
        iDynTree::LinkIndex associatedLinkIndex;
        iDynTree::Transform frameTransform, mixedToBodyTransform;
        std::vector<iDynTree::IndexRange> forcePointsRanges;
        std::vector<iDynTree::Transform> tranformsInFoot;
        std::vector<iDynTree::Wrench> pointForces;
        iDynTree::Wrench totalWrenchInFrame, totalWrenchInLinkBody;
        iDynTree::MatrixDynSize frameJacobianBuffer;
        iDynTree::MatrixFixSize<6,3> pointToLeftJacobianMap;
    } FootVariables;
    FootVariables leftVariables, rightVariables;


    iDynTree::IndexRange basePositionRange, baseQuaternionRange, jointsPositionRange;

    iDynTree::LinkNetExternalWrenches contactWrenches;
    iDynTree::FreeFloatingGeneralizedTorques generalizedStaticTorques;

    iDynTree::VectorDynSize weights, staticTorques;

    double costValue;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;


    void updateVariables (){
        updateRobotState();
        updateFootVariables(leftVariables);
        updateFootVariables(rightVariables);
        computeFootLinkWrench(leftVariables);
        computeFootLinkWrench(rightVariables);
    }


    void prepareFootVariables(const std::string& footName, const iDynTree::FrameIndex & footFrame, const std::vector<iDynTree::Position> &pointsLocalPositions, FootVariables& foot) {
        getFootRanges(footName, foot);
        setFootTransforms(footFrame, pointsLocalPositions, foot);
    }


    void computeFootRelatedJacobians(FootVariables& foot) {

        unsigned int n = static_cast<unsigned int>(jointsPositionRange.size);

        bool ok = sharedKinDyn->getFrameFreeFloatingJacobian(robotState, foot.footFrame, foot.frameJacobianBuffer, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
        assert(ok);

        iDynTree::iDynTreeEigenMatrixMap fullJacobianMap = iDynTree::toEigen(fullJacobianBuffer);
        Eigen::Map<Eigen::Matrix<double, 6, 3, Eigen::RowMajor>> pointToLeftJacMap = iDynTree::toEigen(foot.pointToLeftJacobianMap);
        iDynTree::iDynTreeEigenMatrixMap frameJacobianMap = iDynTree::toEigen(foot.frameJacobianBuffer);
        iDynTree::Transform f_T_a = sharedKinDyn->getWorldTransform(robotState, foot.footFrame).inverse();

        pointToLeftJacMap.topRows<3>() = iDynTree::toEigen(f_T_a.getRotation());

        for (size_t p = 0; p < foot.pointForces.size(); ++p) {
            pointToLeftJacMap.bottomRows<3>() = iDynTree::skew(iDynTree::toEigen(foot.tranformsInFoot[p].getPosition())) * iDynTree::toEigen(f_T_a.getRotation());

            fullJacobianMap.block(0, foot.forcePointsRanges[p].offset, n, 3) = -(frameJacobianMap.rightCols(n)).transpose() * pointToLeftJacMap;
        }

        iDynTree::Transform f_T_b = f_T_a * sharedKinDyn->getBaseTransform(robotState);

        pointToLeftJacMap.topRows<3>() = iDynTree::toEigen(f_T_b.getRotation());

        Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor>> quatDerMap = iDynTree::toEigen(notNormalizedQuaternionMap);

        for (size_t p = 0; p < foot.pointForces.size(); ++p) {
            pointToLeftJacMap.bottomRows<3>() = iDynTree::skew(iDynTree::toEigen(foot.tranformsInFoot[p].getPosition())) * iDynTree::toEigen(f_T_b.getRotation());

            quatDerMap = iDynTree::toEigen(RotatedVectorQuaternionJacobian(foot.pointForces[p].getLinearVec3(), InverseQuaternion(baseQuaternionNormalized))) *
                    iDynTree::toEigen(InverseQuaternionDerivative()) * iDynTree::toEigen(NormalizedQuaternionDerivative(baseQuaternion));

            fullJacobianMap.block(0, baseQuaternionRange.offset, n, 4) -= (frameJacobianMap.rightCols(n)).transpose() * pointToLeftJacMap * quatDerMap;
        }

        pointToLeftJacMap.topRows<3>().setIdentity();

        for (size_t p = 0; p < foot.pointForces.size(); ++p) {
            pointToLeftJacMap.bottomRows<3>() = iDynTree::skew(iDynTree::toEigen(foot.tranformsInFoot[p].getPosition()));

            fullJacobianMap.block(0, jointsPositionRange.offset, n, n) -= (frameJacobianMap.rightCols(n)).transpose() * pointToLeftJacMap *
                    iDynTree::toEigen(RotatedVectorQuaternionJacobian((foot.pointForces[p]).getLinearVec3(),
                                                                      f_T_a.getRotation().asQuaternion())) *
                    iDynTree::toEigen(InverseQuaternionDerivative()) *
                    iDynTree::toEigen(QuaternionLeftTrivializedDerivative(f_T_a.inverse().getRotation().asQuaternion())) *
                    frameJacobianMap.bottomRightCorner(3, n);
        }
    }


private:

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

        iDynTree::toEigen(notNormalizedQuaternionMap) = iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(baseQuaternionNormalized)) *
                iDynTree::toEigen(NormalizedQuaternionDerivative(baseQuaternion));
    }

    void updateFootVariables(FootVariables& foot) {
        iDynTree::Vector3 pointForce;
        for (size_t p = 0; p < foot.forcePointsRanges.size(); ++p) {
            pointForce = stateVariables(foot.forcePointsRanges[p]);
            foot.pointForces[p].setLinearVec3(pointForce);
        }
    }

    void getFootRanges(const std::string& footName, FootVariables& foot) {

        size_t forcePoints = 0;
        for (auto& label : stateVariables.listOfLabels()) {
            if (label.find(footName + "ForcePoint") != std::string::npos) {
                ++forcePoints;
            }
        }

        foot.forcePointsRanges.resize(forcePoints);
        foot.pointForces.resize(forcePoints);

        for (size_t p = 0; p < forcePoints; ++p) {
            foot.forcePointsRanges[p] = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(p));
            assert(foot.forcePointsRanges[p].isValid());
            foot.pointForces[p].zero();
        }

    }

    void setFootTransforms(const iDynTree::FrameIndex & footFrame, const std::vector<iDynTree::Position> &pointsLocalPositions, FootVariables& foot) {
        assert(pointsLocalPositions.size() == foot.forcePointsRanges.size());

        foot.tranformsInFoot.clear();
        for (auto& position : pointsLocalPositions) {
            foot.tranformsInFoot.push_back(iDynTree::Transform(iDynTree::Rotation::Identity(), position));
        }

        assert(footFrame != iDynTree::FRAME_INVALID_INDEX);
        assert(timedSharedKinDyn->model().isValidFrameIndex(footFrame));

        foot.footFrame = footFrame;
        foot.associatedLinkIndex = timedSharedKinDyn->model().getFrameLink(footFrame);
        foot.frameTransform = timedSharedKinDyn->model().getFrameTransform(footFrame);
        foot.mixedToBodyTransform = iDynTree::Transform::Identity();

        foot.frameJacobianBuffer.resize(6, static_cast<unsigned int>(6 + jointsPositionRange.size));
        foot.frameJacobianBuffer.zero();
    }

    void computeFootLinkWrench(FootVariables& foot) {
        foot.totalWrenchInFrame.zero();

        foot.mixedToBodyTransform.setRotation(sharedKinDyn->getWorldTransform(robotState, foot.footFrame).getRotation().inverse());

        for (size_t p = 0; p < foot.pointForces.size(); ++p) {
            foot.totalWrenchInFrame = foot.totalWrenchInFrame + foot.tranformsInFoot[p] * (foot.mixedToBodyTransform * foot.pointForces[p]);
        }

        contactWrenches(foot.associatedLinkIndex) = foot.frameTransform * foot.totalWrenchInFrame;
    }

};



StaticTorquesCost::StaticTorquesCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn,
                                     const iDynTree::FrameIndex &leftFootFrame, const iDynTree::FrameIndex &rightFootFrame, const std::vector<iDynTree::Position> &positionsInLeftFoot,
                                     const std::vector<iDynTree::Position> &positionsInRightFoot)
    : iDynTree::optimalcontrol::Cost ("StaticTorques")
    , m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());
    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->prepareFootVariables("Left", leftFootFrame, positionsInLeftFoot, m_pimpl->leftVariables);
    m_pimpl->prepareFootVariables("Right", rightFootFrame, positionsInRightFoot, m_pimpl->rightVariables);

    m_pimpl->contactWrenches.resize(timelySharedKinDyn->model());
    m_pimpl->generalizedStaticTorques.resize(timelySharedKinDyn->model());

    unsigned int n = static_cast<unsigned int>(m_pimpl->jointsPositionRange.size);
    m_pimpl->weights.resize(n);
    iDynTree::toEigen(m_pimpl->weights).setConstant(1.0);
    m_pimpl->staticTorques.resize(n);

    m_pimpl->stateGradientBuffer.resize(static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->stateGradientBuffer.zero();
    m_pimpl->controlGradientBuffer.resize(static_cast<unsigned int>(controlVariables.size()));
    m_pimpl->controlGradientBuffer.zero();
    m_pimpl->fullJacobianBuffer.resize(n, static_cast<unsigned int>(stateVariables.size()));
    m_pimpl->fullJacobianBuffer.zero();
    m_pimpl->jointsJacobianBuffer.resize(n,n);
    m_pimpl->jointsJacobianBuffer.zero();
    m_pimpl->massMatrixBuffer.resize(n+6, n+6);
    m_pimpl->massMatrixBuffer.zero();

}

StaticTorquesCost::~StaticTorquesCost()
{ }

bool StaticTorquesCost::setWeights(const iDynTree::VectorDynSize &torquesWeights)
{
    if (torquesWeights.size() != m_pimpl->weights.size()) {
        return false;
    }
    m_pimpl->weights = torquesWeights;

    return true;
}

void StaticTorquesCost::computeStaticTorques(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &staticTorques)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateVariables();

    bool ok = m_pimpl->sharedKinDyn->getStaticForces(m_pimpl->robotState, m_pimpl->contactWrenches, m_pimpl->generalizedStaticTorques);
    assert(ok);

    staticTorques = m_pimpl->generalizedStaticTorques.jointTorques();
}

bool StaticTorquesCost::costEvaluation(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, double &costValue)
{

    computeStaticTorques(time, state, control, m_pimpl->staticTorques);

    m_pimpl->costValue = 0.5 * iDynTree::toEigen(m_pimpl->staticTorques).transpose() * iDynTree::toEigen(m_pimpl->weights).asDiagonal() * iDynTree::toEigen(m_pimpl->staticTorques);

    costValue = m_pimpl->costValue;

    return true;
}

void StaticTorquesCost::computeStaticTorquesJacobian(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &staticTorquesJacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateVariables();

    unsigned int n = static_cast<unsigned int>(m_pimpl->jointsPositionRange.size);

    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->fullJacobianBuffer);

    bool ok = m_pimpl->sharedKinDyn->getStaticForcesJointsDerivative(m_pimpl->robotState, m_pimpl->contactWrenches, m_pimpl->jointsJacobianBuffer);
    assert(ok);

    jacobianMap.block(0, m_pimpl->jointsPositionRange.offset, n, n) = iDynTree::toEigen(m_pimpl->jointsJacobianBuffer);


    ok = m_pimpl->sharedKinDyn->getFreeFloatingMassMatrix(m_pimpl->robotState, m_pimpl->massMatrixBuffer, iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);
    assert(ok);

    jacobianMap.block(0, m_pimpl->baseQuaternionRange.offset, n, 4) = -iDynTree::toEigen(m_pimpl->massMatrixBuffer).bottomLeftCorner(n, 3) *
            iDynTree::toEigen(RotatedVectorQuaternionJacobian(m_pimpl->sharedKinDyn->gravity(), InverseQuaternion(m_pimpl->baseQuaternionNormalized))) *
            iDynTree::toEigen(InverseQuaternionDerivative()) * iDynTree::toEigen(NormalizedQuaternionDerivative(m_pimpl->baseQuaternion));

    m_pimpl->computeFootRelatedJacobians(m_pimpl->leftVariables);

    m_pimpl->computeFootRelatedJacobians(m_pimpl->rightVariables);

    staticTorquesJacobian = m_pimpl->fullJacobianBuffer;
}

bool StaticTorquesCost::costFirstPartialDerivativeWRTState(double time, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{
    unsigned int n = static_cast<unsigned int>(m_pimpl->jointsPositionRange.size);

    computeStaticTorques(time, state, control, m_pimpl->staticTorques);
    computeStaticTorquesJacobian(time, state, control, m_pimpl->fullJacobianBuffer);

    iDynTree::iDynTreeEigenVector gradientMap = iDynTree::toEigen(m_pimpl->stateGradientBuffer);
    iDynTree::iDynTreeEigenMatrixMap fullJacobianMap = iDynTree::toEigen(m_pimpl->fullJacobianBuffer);

    gradientMap.segment(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsPositionRange.size) = (iDynTree::toEigen(m_pimpl->staticTorques).transpose() *
                                                                                                   iDynTree::toEigen(m_pimpl->weights).asDiagonal() *
                                                                                                   fullJacobianMap.block(0, m_pimpl->jointsPositionRange.offset, n, n)).transpose();

    gradientMap.segment<4>(m_pimpl->baseQuaternionRange.offset) = (iDynTree::toEigen(m_pimpl->staticTorques).transpose() *
                                                                   iDynTree::toEigen(m_pimpl->weights).asDiagonal() *
                                                                   fullJacobianMap.block(0, m_pimpl->baseQuaternionRange.offset, n, 4)).transpose();


    for (auto leftRange: m_pimpl->leftVariables.forcePointsRanges) {
        gradientMap.segment<3>(leftRange.offset) = (iDynTree::toEigen(m_pimpl->staticTorques).transpose() *
                                                    iDynTree::toEigen(m_pimpl->weights).asDiagonal() *
                                                    fullJacobianMap.block(0, leftRange.offset, n, 3)).transpose();
    }

    for (auto rightRange: m_pimpl->rightVariables.forcePointsRanges) {
        gradientMap.segment<3>(rightRange.offset) = (iDynTree::toEigen(m_pimpl->staticTorques).transpose() *
                                                     iDynTree::toEigen(m_pimpl->weights).asDiagonal() *
                                                     fullJacobianMap.block(0, rightRange.offset, n, 3)).transpose();
    }

    partialDerivative = m_pimpl->stateGradientBuffer;
    return true;
}

bool StaticTorquesCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{
    partialDerivative = m_pimpl->controlGradientBuffer;
    return true;
}
