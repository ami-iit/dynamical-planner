/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

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
    iDynTree::MatrixDynSize fullJacobianBuffer, jointsJacobianBuffer;

    typedef struct {
        iDynTree::FrameIndex footFrame;
        iDynTree::LinkIndex associatedLinkIndex;
        iDynTree::Transform frameTransform, mixedToBodyTransform;
        std::vector<iDynTree::IndexRange> forcePointsRanges;
        std::vector<iDynTree::Transform> tranformsInFoot;
        std::vector<iDynTree::Wrench> pointForces;
        iDynTree::Wrench totalWrenchInFrameMixed, totalWrenchInLinkBody;
    } FootVariables;
    FootVariables leftVariables, rightVariables;

    typedef struct {
        iDynTree::Matrix3x3 gravityForceSkew;
        iDynTree::Position comPosition;
        iDynTree::MatrixFixSize<3, 6> transformationToCoM;
        iDynTree::MatrixDynSize jacobianBuffer, forceDerivative;
        iDynTree::SpatialMotionVector parentJointMotionSubspaceVector;
    } LinkVariables;

    std::vector<LinkVariables> linksVariables;

    iDynTree::IndexRange basePositionRange, baseQuaternionRange, jointsPositionRange;

    iDynTree::LinkNetExternalWrenches contactWrenches;
    iDynTree::FreeFloatingGeneralizedTorques generalizedStaticTorques;
    iDynTree::LinkWrenches linkStaticForces;

    iDynTree::VectorDynSize weights, staticTorques;

    double costValue;

    RobotState robotState;
    std::shared_ptr<SharedKinDynComputation> sharedKinDyn;


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

    void prepareLinksVariables() {

        linksVariables.resize(sharedKinDyn->model().getNrOfLinks());

        Eigen::Map<const Eigen::Vector3d> gravityMap = iDynTree::toEigen(sharedKinDyn->gravity());

        for (size_t l = 0; l < sharedKinDyn->model().getNrOfLinks(); ++l) {
            iDynTree::LinkIndex linkIndex = static_cast<iDynTree::LinkIndex>(l);
            iDynTree::LinkConstPtr linkPtr = sharedKinDyn->model().getLink(linkIndex);
            iDynTree::toEigen(linksVariables[l].gravityForceSkew) = -iDynTree::skew(linkPtr->getInertia().getMass() * gravityMap);
            linksVariables[l].comPosition = linkPtr->getInertia().getCenterOfMass();
            linksVariables[l].jacobianBuffer.resize(6, 6 + static_cast<unsigned int>(jointsPositionRange.size));
            linksVariables[l].forceDerivative.resize(3, 6 + static_cast<unsigned int>(jointsPositionRange.size));
            linksVariables[l].transformationToCoM.zero();
            iDynTree::toEigen(linksVariables[l].transformationToCoM).leftCols<3>().setIdentity();

            const iDynTree::Link* parentLinkPtr = sharedKinDyn->traversal().getParentLinkFromLinkIndex(linkIndex);

            if (parentLinkPtr != nullptr) {
                iDynTree::LinkIndex parentLink = sharedKinDyn->traversal().getParentLinkFromLinkIndex(linkIndex)->getIndex();
                linksVariables[l].parentJointMotionSubspaceVector = sharedKinDyn->traversal().getParentJointFromLinkIndex(linkIndex)->getMotionSubspaceVector(0, linkIndex, parentLink);
            }
        }
    }

    void updateLinksJacobians() {
        bool ok = false;
        iDynTree::Transform linkTransform;
        for (size_t l = 0; l < sharedKinDyn->model().getNrOfLinks(); ++l) {
            iDynTree::LinkIndex linkIndex = static_cast<iDynTree::LinkIndex>(l);

            ok = sharedKinDyn->getFrameFreeFloatingJacobian(robotState, linkIndex, linksVariables[l].jacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);
            assert(ok);

            linkTransform = sharedKinDyn->getWorldTransform(robotState, linkIndex);

            iDynTree::toEigen(linksVariables[l].transformationToCoM).rightCols<3>() = iDynTree::skew(iDynTree::toEigen(linkTransform.getPosition() - (linkTransform * linksVariables[l].comPosition)));

            iDynTree::toEigen(linksVariables[l].forceDerivative) = iDynTree::toEigen(linksVariables[l].gravityForceSkew) * iDynTree::toEigen(linksVariables[l].transformationToCoM) *
                    iDynTree::toEigen(linksVariables[l].jacobianBuffer);

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

        robotState.world_T_base.setRotation(baseRotation);
        robotState.world_T_base.setPosition(basePosition);

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
        for (auto label : stateVariables.listOfLabels()) {
            if (label.find(footName + "ForcePoint") != std::string::npos) {
                ++forcePoints;
            }
        }

        foot.forcePointsRanges.clear();
        foot.pointForces.resize(forcePoints);

        for (size_t p = 0; p < forcePoints; ++p) {
            foot.forcePointsRanges.push_back(stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(p)));
            assert(foot.forcePointsRanges.back().isValid());
            foot.pointForces[p].zero();
        }

    }

    void setFootTransforms(const iDynTree::FrameIndex & footFrame, const std::vector<iDynTree::Position> &pointsLocalPositions, FootVariables& foot) {
        assert(pointsLocalPositions.size() == foot.forcePointsRanges.size());

        foot.tranformsInFoot.clear();
        for (auto position : pointsLocalPositions) {
            foot.tranformsInFoot.push_back(iDynTree::Transform(iDynTree::Rotation::Identity(), position));
        }

        assert(footFrame != iDynTree::FRAME_INVALID_INDEX);
        assert(sharedKinDyn->model().isValidFrameIndex(footFrame));

        foot.footFrame = footFrame;
        foot.associatedLinkIndex = sharedKinDyn->model().getFrameLink(footFrame);
        foot.frameTransform = sharedKinDyn->model().getFrameTransform(footFrame);
        foot.mixedToBodyTransform = iDynTree::Transform::Identity();
    }

    void computeFootLinkWrench(FootVariables& foot) {
        foot.totalWrenchInFrameMixed.zero();

        for (size_t p = 0; p < foot.pointForces.size(); ++p) {
            foot.totalWrenchInFrameMixed = foot.totalWrenchInFrameMixed + foot.tranformsInFoot[p] * foot.pointForces[p];
        }

        foot.mixedToBodyTransform.setRotation(sharedKinDyn->getWorldTransform(robotState, foot.footFrame).getRotation().inverse());

        contactWrenches(foot.associatedLinkIndex) = foot.frameTransform * foot.mixedToBodyTransform * foot.totalWrenchInFrameMixed;
    }

};



StaticTorquesCost::StaticTorquesCost(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, std::shared_ptr<SharedKinDynComputation> sharedKinDyn,
                                     const iDynTree::FrameIndex &leftFootFrame, const iDynTree::FrameIndex &rightFootFrame, const std::vector<iDynTree::Position> &positionsInLeftFoot,
                                     const std::vector<iDynTree::Position> &positionsInRightFoot)
    : iDynTree::optimalcontrol::Cost ("StaticTorques")
    , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;

    assert(sharedKinDyn);
    assert(sharedKinDyn->isValid());
    m_pimpl->sharedKinDyn = sharedKinDyn;

    m_pimpl->prepareFootVariables("Left", leftFootFrame, positionsInLeftFoot, m_pimpl->leftVariables);
    m_pimpl->prepareFootVariables("Right", rightFootFrame, positionsInRightFoot, m_pimpl->rightVariables);
    m_pimpl->prepareLinksVariables();

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

    m_pimpl->robotState = sharedKinDyn->currentState();
    m_pimpl->contactWrenches.resize(sharedKinDyn->model());
    m_pimpl->generalizedStaticTorques.resize(sharedKinDyn->model());
    m_pimpl->linkStaticForces.resize(sharedKinDyn->model());

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

void StaticTorquesCost::computeStaticTorques(const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &staticTorques)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->updateVariables();

    bool ok = m_pimpl->sharedKinDyn->getStaticForces(m_pimpl->robotState, m_pimpl->contactWrenches, m_pimpl->generalizedStaticTorques, m_pimpl->linkStaticForces);
    assert(ok);

    staticTorques = m_pimpl->generalizedStaticTorques.jointTorques();
}

bool StaticTorquesCost::costEvaluation(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, double &costValue)
{

    computeStaticTorques(state, control, m_pimpl->staticTorques);

    m_pimpl->costValue = 0.5 * iDynTree::toEigen(m_pimpl->staticTorques).transpose() * iDynTree::toEigen(m_pimpl->weights).asDiagonal() * iDynTree::toEigen(m_pimpl->staticTorques);

    costValue = m_pimpl->costValue;

    return true;
}

void StaticTorquesCost::computeStaticTorquesJacobian(const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::MatrixDynSize &staticTorquesJacobian)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = control;

    m_pimpl->updateVariables();

    m_pimpl->updateLinksJacobians();



}

bool StaticTorquesCost::costFirstPartialDerivativeWRTState(double, const iDynTree::VectorDynSize &state, const iDynTree::VectorDynSize &control, iDynTree::VectorDynSize &partialDerivative)
{

}

bool StaticTorquesCost::costFirstPartialDerivativeWRTControl(double, const iDynTree::VectorDynSize &, const iDynTree::VectorDynSize &, iDynTree::VectorDynSize &partialDerivative)
{

}
