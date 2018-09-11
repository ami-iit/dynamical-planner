/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
#include <iDynTree/Model/ForwardKinematics.h>
#include <iDynTree/Model/Dynamics.h>
#include <DynamicalPlannerPrivate/CheckEqualVector.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

bool SharedKinDynComputation::sameState(const RobotState &other)
{
    if (VectorsAreEqual(other.world_T_base.getPosition(), m_state.world_T_base.getPosition(), m_tol)
            && VectorsAreEqual(other.world_T_base.getRotation().asQuaternion(), m_state.world_T_base.getRotation().asQuaternion(), m_tol)
            && VectorsAreEqual(other.s, m_state.s, m_tol)
            && VectorsAreEqual(other.base_velocity, m_state.base_velocity, m_tol)
            && VectorsAreEqual(other.s_dot, m_state.s_dot, m_tol)) {
        return true;
    }
    return false;
}

bool SharedKinDynComputation::updateRobotState(const RobotState &currentState)
{
    if (m_updateNecessary || !sameState(currentState)) {
        m_kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION); //The base_velocity saved in the robot state is supposed to be in body frame
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        if (!ok) {
            return false;
        }
        m_updateNecessary = false;
        m_state = currentState;
    }
    return true;
}

void SharedKinDynComputation::fillJointsInfo()
{
    const iDynTree::Model& model = m_kinDyn.model();
    m_childrenForceDerivatives.resize(model.getNrOfLinks());
    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {
        iDynTree::JointIndex jointIndex = static_cast<iDynTree::JointIndex>(j);
        assert(model.isValidJointIndex(jointIndex));
        m_jointsInfos[j].jointPtr = model.getJoint(jointIndex);
        assert(m_jointsInfos[j].jointPtr->getNrOfDOFs() == 1);
        m_jointsInfos[j].childIndex =  m_traversal.getChildLinkIndexFromJointIndex(model, jointIndex);
        m_jointsInfos[j].parentIndex =  m_traversal.getParentLinkIndexFromJointIndex(model, jointIndex);
        size_t childIndex = static_cast<size_t>(m_jointsInfos[j].childIndex);
        m_childrenForceDerivatives[childIndex].resize(m_jointsInfos.size(), iDynTree::SpatialForceVector::Zero());
    }
}

void SharedKinDynComputation::updateChildBuffersForMomentumDerivative()
{
    iDynTree::LinkIndex baseIndex = m_kinDyn.model().getLinkIndex(m_kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {

        m_jointsInfos[j].motionVectorTimesChildVelocity = m_jointsInfos[j].jointPtr->getMotionSubspaceVector(0,
                                                                                                m_jointsInfos[j].childIndex,
                                                                                                m_jointsInfos[j].parentIndex).cross(
                    -m_kinDyn.getFrameVel(m_jointsInfos[j].childIndex));

        m_jointsInfos[j].childVelocity = m_kinDyn.getFrameVel(m_jointsInfos[j].childIndex);

        m_jointsInfos[j].baseTC = m_kinDyn.getRelativeTransform(baseIndex, m_jointsInfos[j].childIndex);
    }
}

void SharedKinDynComputation::computeChildStaticForceDerivative(const iDynTree::LinkWrenches &linkStaticForces)
{

    iDynTree::LinkIndex childLink, parentLink;
    size_t childIndex;
    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {

        childLink = m_jointsInfos[j].childIndex;
        parentLink = m_jointsInfos[j].parentIndex;

        JointInfos& jInfo = m_jointsInfos[j];

        jInfo.childStaticForceDerivative = jInfo.jointPtr->getMotionSubspaceVector(0,
                                                                                   childLink,
                                                                                   parentLink).cross(
                    linkStaticForces(childLink));

        jInfo.motionVectorTimesChildAcceleration = jInfo.jointPtr->getMotionSubspaceVector(0,
                                                                                           childLink,
                                                                                           parentLink).cross(
                    m_invDynLinkProperAccs(childLink));

        childIndex = static_cast<size_t>(childLink);
        m_childrenForceDerivatives[childIndex] = m_zeroDerivatives;
    }
}

bool SharedKinDynComputation::computeStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces)
{
    iDynTree::toEigen(m_gravityAccInBaseLinkFrame) =
            iDynTree::toEigen(currentState.world_T_base.getRotation().inverse())*iDynTree::toEigen(m_gravity);

    // Clear input buffers that need to be cleared
    m_invDynGeneralizedProperAccs.baseAcc().zero();
    iDynTree::toEigen(m_invDynGeneralizedProperAccs.baseAcc().getLinearVec3()) = - iDynTree::toEigen(m_gravityAccInBaseLinkFrame);
    m_invDynGeneralizedProperAccs.jointAcc().zero();

    m_pos.worldBasePos() = currentState.world_T_base;
    iDynTree::toEigen(m_pos.jointPos()) = iDynTree::toEigen(currentState.s);

    // Run inverse dynamics
    bool ok = iDynTree::ForwardAccKinematics(m_kinDyn.model(),
                                             m_traversal,
                                             m_pos,
                                             m_invDynZeroVel,
                                             m_invDynGeneralizedProperAccs,
                                             m_invDynZeroLinkVel,
                                             m_invDynLinkProperAccs);

    ok = ok && iDynTree::RNEADynamicPhase(m_kinDyn.model(),
                                          m_traversal,
                                          m_pos.jointPos(),
                                          m_invDynZeroLinkVel,
                                          m_invDynLinkProperAccs,
                                          linkExtForces,
                                          m_linkStaticWrenches,
                                          m_generalizedStaticTorques);
    return ok;
}

SharedKinDynComputation::SharedKinDynComputation()
{
    m_state.world_T_base.setRotation(iDynTree::Rotation::Identity());
    m_state.world_T_base.setPosition(iDynTree::Position::Zero());
    m_state.base_velocity.zero();
    m_gravity.zero();
    m_gravity(2) = -9.81;
    m_updateNecessary = true;
    m_tol = iDynTree::DEFAULT_TOL;
}

bool SharedKinDynComputation::loadRobotModel(const iDynTree::Model &model)
{
    std::lock_guard<std::mutex> guard(m_mutex);
    bool ok = m_kinDyn.loadRobotModel(model);

    if (ok) {
        m_state.s.resize(m_kinDyn.getNrOfDegreesOfFreedom());
        m_state.s.zero();
        m_state.s_dot.resize(m_kinDyn.getNrOfDegreesOfFreedom());
        m_state.s_dot.zero();
        m_jointsInfos.resize(m_kinDyn.getNrOfDegreesOfFreedom());
        bool okTraversal = m_kinDyn.model().computeFullTreeTraversal(m_traversal);
        if (!okTraversal)
            return false;
        fillJointsInfo();
    }

    m_linkStaticWrenches.resize(model);
    m_invDynGeneralizedProperAccs.resize(model);
    m_pos.resize(model);
    m_invDynZeroVel.resize(model);
    m_invDynZeroVel.baseVel().zero();
    m_invDynZeroVel.jointVel().zero();
    m_invDynZeroLinkVel.resize(model);
    m_invDynLinkProperAccs.resize(model);

    for(iDynTree::LinkIndex lnkIdx = 0; lnkIdx < static_cast<iDynTree::LinkIndex>(model.getNrOfLinks()); lnkIdx++)
    {
        m_invDynZeroLinkVel(lnkIdx).zero();
    }

    m_generalizedStaticTorques.resize(model);
    m_zeroDerivatives.resize(m_kinDyn.getNrOfDegreesOfFreedom(), iDynTree::SpatialForceVector::Zero());

    return ok;
}

const iDynTree::Model &SharedKinDynComputation::model() const
{
    return m_kinDyn.model();
}

bool SharedKinDynComputation::isValid() const
{
    return m_kinDyn.isValid();
}

void SharedKinDynComputation::setGravity(const iDynTree::Vector3 &gravity)
{
    std::lock_guard<std::mutex> guard(m_mutex);
    if (!VectorsAreEqual(gravity, m_gravity, m_tol))
        m_updateNecessary = true;

    m_gravity = gravity;
}

const iDynTree::Vector3 &SharedKinDynComputation::gravity() const
{
    return m_gravity;
}

bool SharedKinDynComputation::setToleranceForUpdate(double tol)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (tol < 0) {
        return false;
    }

    m_tol = tol;

    return true;
}

double SharedKinDynComputation::getUpdateTolerance() const
{
    return m_tol;
}

bool SharedKinDynComputation::setFloatingBase(const std::string &floatingBaseName)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_updateNecessary = true;

    if (!m_kinDyn.setFloatingBase(floatingBaseName)) {
        return false;
    }

    bool okTraversal = m_kinDyn.model().computeFullTreeTraversal(m_traversal, m_kinDyn.model().getLinkIndex(m_kinDyn.getFloatingBase()));
    if (!okTraversal)
        return false;

    fillJointsInfo();

    return true;
}

const iDynTree::Traversal &SharedKinDynComputation::traversal() const
{
    return m_traversal;
}

const RobotState &SharedKinDynComputation::currentState() const
{
    return m_state;
}

iDynTree::Position SharedKinDynComputation::getCenterOfMassPosition(const RobotState &currentState)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getCenterOfMassPosition();
}

bool SharedKinDynComputation::getCenterOfMassJacobian(const RobotState &currentState, iDynTree::MatrixDynSize &comJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (!updateRobotState(currentState))
        return false;

    return m_kinDyn.getCenterOfMassJacobian(comJacobian);

}

iDynTree::Transform SharedKinDynComputation::getWorldTransform(const RobotState &currentState, std::string frameName)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getWorldTransform(frameName);
}

iDynTree::Transform SharedKinDynComputation::getWorldTransform(const RobotState &currentState, const iDynTree::FrameIndex frameIndex)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getWorldTransform(frameIndex);
}

bool SharedKinDynComputation::getFrameFreeFloatingJacobian(const RobotState &currentState, const std::string &frameName, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (!updateRobotState(currentState))
        return false;

    return m_kinDyn.getFrameFreeFloatingJacobian(frameName, outJacobian);
}

bool SharedKinDynComputation::getFrameFreeFloatingJacobian(const RobotState &currentState, const iDynTree::FrameIndex frameIndex, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (!updateRobotState(currentState))
        return false;

    return m_kinDyn.getFrameFreeFloatingJacobian(frameIndex, outJacobian);
}

iDynTree::Transform SharedKinDynComputation::getRelativeTransform(const RobotState &currentState,
                                                                  const iDynTree::FrameIndex refFrameIndex,
                                                                  const iDynTree::FrameIndex frameIndex)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getRelativeTransform(refFrameIndex, frameIndex);
}

iDynTree::Transform SharedKinDynComputation::getRelativeTransform(const RobotState &currentState,
                                                                  const std::string &refFrameName,
                                                                  const std::string &frameName)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getRelativeTransform(refFrameName, frameName);
}

bool SharedKinDynComputation::getRelativeJacobian(const RobotState &currentState, const iDynTree::FrameIndex refFrameIndex, const iDynTree::FrameIndex frameIndex, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (!updateRobotState(currentState))
        return false;

    return m_kinDyn.getRelativeJacobian(refFrameIndex, frameIndex, outJacobian);
}

iDynTree::Twist SharedKinDynComputation::getFrameVel(const RobotState &currentState, const std::string &frameName, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getFrameVel(frameName);
}

iDynTree::Twist SharedKinDynComputation::getFrameVel(const RobotState &currentState, const iDynTree::FrameIndex frameIdx, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getFrameVel(frameIdx);
}

bool SharedKinDynComputation::getFrameVelJointsDerivative(const RobotState &currentState, const iDynTree::FrameIndex frameIdx, iDynTree::MatrixDynSize &velocityDerivative)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!m_kinDyn.isValid())
        return false;

    m_kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    if (!updateRobotState(currentState))
        return false;

    velocityDerivative.resize(6, static_cast<unsigned int>(m_jointsInfos.size()));
    velocityDerivative.zero();
    iDynTree::iDynTreeEigenMatrixMap derivativeMap = iDynTree::toEigen(velocityDerivative);

    const iDynTree::Model& model = m_kinDyn.model();
    iDynTree::LinkIndex linkIndex = model.getFrameLink(frameIdx);
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    if (linkIndex == iDynTree::LINK_INVALID_INDEX) {
        return false;
    }

    iDynTree::IJointConstPtr jointPtr;
    size_t jointIndex;

    iDynTree::LinkIndex visitedLink = linkIndex;
    iDynTree::Twist jointDerivative, childVelocity;
    while (visitedLink != baseIndex) {
        jointPtr = m_traversal.getParentJointFromLinkIndex(visitedLink);
        jointIndex = static_cast<size_t>(jointPtr->getIndex());

        childVelocity = m_kinDyn.getFrameVel(m_jointsInfos[jointIndex].childIndex);
        jointDerivative = m_kinDyn.getRelativeTransform(frameIdx, m_jointsInfos[jointIndex].childIndex) *
                m_jointsInfos[jointIndex].jointPtr->getMotionSubspaceVector(0,
                                                                            m_jointsInfos[jointIndex].childIndex,
                                                                            m_jointsInfos[jointIndex].parentIndex).cross(-childVelocity);

        derivativeMap.col(static_cast<Eigen::Index>(jointIndex)) = iDynTree::toEigen(jointDerivative);
        visitedLink = m_traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
    }

    return true;
}

iDynTree::SpatialMomentum SharedKinDynComputation::getLinearAngularMomentum(const RobotState &currentState, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_kinDyn.getLinearAngularMomentum();
}

bool SharedKinDynComputation::getLinearAngularMomentumJacobian(const RobotState &currentState, iDynTree::MatrixDynSize &linAngMomentumJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (!updateRobotState(currentState))
        return false;

    return m_kinDyn.getLinearAngularMomentumJacobian(linAngMomentumJacobian);
}

bool SharedKinDynComputation::getLinearAngularMomentumJointsDerivative(const RobotState &currentState,
                                                                       iDynTree::MatrixDynSize &linAngMomentumDerivative)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!m_kinDyn.isValid())
        return false;

    m_kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    if (!updateRobotState(currentState))
        return false;

    linAngMomentumDerivative.resize(6, static_cast<unsigned int>(m_jointsInfos.size()));
    linAngMomentumDerivative.zero();

    updateChildBuffersForMomentumDerivative();

    iDynTree::iDynTreeEigenMatrixMap derivativeMap = iDynTree::toEigen(linAngMomentumDerivative);

    const iDynTree::Model& model = m_kinDyn.model();
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    iDynTree::IJointConstPtr jointPtr;
    iDynTree::LinkIndex linkIndex;
    iDynTree::LinkConstPtr linkPtr;
    size_t jointIndex;
    iDynTree::SpatialMomentum linkMomentum, linkMomentumInChildFrame, jointMomentumDerivative;
    iDynTree::SpatialMomentum transformDerivative, velocityDerivative;
    iDynTree::LinkIndex visitedLink, childLink, parentLink;
    iDynTree::Transform b_T_link;
    iDynTree::Twist childVelocity;

    for (size_t l = 0; l < model.getNrOfLinks(); ++l) {
        linkIndex = static_cast<iDynTree::LinkIndex>(l);
        assert(model.isValidLinkIndex(linkIndex));
        linkPtr = model.getLink(linkIndex);
        linkMomentum = linkPtr->getInertia() * m_kinDyn.getFrameVel(linkIndex);
        b_T_link = m_kinDyn.getRelativeTransform(baseIndex, linkIndex);

        visitedLink = linkIndex;
        while(visitedLink != baseIndex) {
            jointPtr = m_traversal.getParentJointFromLinkIndex(visitedLink);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());

            childLink = m_jointsInfos[jointIndex].childIndex;
            parentLink = m_jointsInfos[jointIndex].parentIndex;

            childVelocity = m_jointsInfos[jointIndex].childVelocity;
            linkMomentumInChildFrame = m_kinDyn.getRelativeTransform(childLink, linkIndex) * linkMomentum;

            transformDerivative = m_jointsInfos[jointIndex].baseTC *
                    jointPtr->getMotionSubspaceVector(0,
                                                      childLink,
                                                      parentLink).cross(linkMomentumInChildFrame);

            velocityDerivative = b_T_link * (linkPtr->getInertia() * (m_kinDyn.getRelativeTransform(linkIndex, childLink) *
                                                                      m_jointsInfos[jointIndex].motionVectorTimesChildVelocity));

            jointMomentumDerivative = transformDerivative + velocityDerivative;

            derivativeMap.col(static_cast<Eigen::Index>(jointIndex)) += iDynTree::toEigen(jointMomentumDerivative);

            visitedLink = m_traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
        }
    }
    return true;
}

bool SharedKinDynComputation::getStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces, iDynTree::FreeFloatingGeneralizedTorques& generalizedStaticTorques, iDynTree::LinkWrenches& linkStaticForces)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!computeStaticForces(currentState, linkExtForces)) {
        return false;
    }

    linkStaticForces = m_linkStaticWrenches;
    generalizedStaticTorques = m_generalizedStaticTorques;

    return true;
}

bool SharedKinDynComputation::getStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces, iDynTree::FreeFloatingGeneralizedTorques &generalizedStaticTorques)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!computeStaticForces(currentState, linkExtForces)) {
        return false;
    }

    generalizedStaticTorques = m_generalizedStaticTorques;

    return true;
}


bool SharedKinDynComputation::getStaticForcesJointsDerivative(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces, iDynTree::MatrixDynSize &staticTorquesDerivatives)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!m_kinDyn.isValid())
        return false;

    m_kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    if (!updateRobotState(currentState))
        return false;

    if (!computeStaticForces(currentState, linkExtForces)) {
        return false;
    }

    computeChildStaticForceDerivative(m_linkStaticWrenches);

    staticTorquesDerivatives.resize(static_cast<unsigned int>(m_jointsInfos.size()), static_cast<unsigned int>(m_jointsInfos.size()));
    staticTorquesDerivatives.zero();

    const iDynTree::Model& model = m_kinDyn.model();
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    iDynTree::LinkIndex parentLinkIndex, associatedLinkIndex;
    iDynTree::LinkConstPtr ancestor_ptr;
    iDynTree::TraversalIndex elIndex;
    iDynTree::IJointConstPtr visitedJoint, associatedJoint_ptr;
    iDynTree::SpatialInertia linkInertia;
    size_t visitedJointIndex, associatedJoint, parentLink, associatedLink, ancestorLink;
    iDynTree::Transform l_T_c, p_T_c, ancestor_T_associated;

    for (unsigned int el = m_traversal.getNrOfVisitedLinks() -1; el > 0; el--) {

        elIndex = static_cast<iDynTree::TraversalIndex>(el);
        associatedJoint_ptr = m_traversal.getParentJoint(elIndex);
        associatedJoint = static_cast<size_t>(associatedJoint_ptr->getIndex());
        associatedLinkIndex = m_traversal.getLink(elIndex)->getIndex();
        linkInertia = model.getLink(associatedLinkIndex)->getInertia();
        associatedLink = static_cast<size_t>(associatedLinkIndex);

        parentLinkIndex = m_traversal.getParentLink(elIndex)->getIndex();
        parentLink = static_cast<size_t>(parentLinkIndex);

        p_T_c = m_kinDyn.getRelativeTransform(parentLinkIndex, associatedLinkIndex);

        m_childrenForceDerivatives[associatedLink][associatedJoint] = m_childrenForceDerivatives[associatedLink][associatedJoint] + m_jointsInfos[associatedJoint].childStaticForceDerivative;

        // Propagate joints which are before the joint
        visitedJoint = associatedJoint_ptr;
        while (visitedJoint) {
            visitedJointIndex = static_cast<size_t>(visitedJoint->getIndex());

            l_T_c = m_kinDyn.getRelativeTransform(associatedLinkIndex, m_jointsInfos[visitedJointIndex].childIndex);
            m_childrenForceDerivatives[associatedLink][visitedJointIndex] = m_childrenForceDerivatives[associatedLink][visitedJointIndex] - linkInertia * (l_T_c * m_jointsInfos[visitedJointIndex].motionVectorTimesChildAcceleration);

            if (parentLinkIndex != baseIndex) {
                m_childrenForceDerivatives[parentLink][visitedJointIndex] = m_childrenForceDerivatives[parentLink][visitedJointIndex] + (p_T_c * m_childrenForceDerivatives[associatedLink][visitedJointIndex]);
            }
            visitedJoint = m_traversal.getParentJointFromLinkIndex(m_jointsInfos[visitedJointIndex].parentIndex);
        }

        if (parentLinkIndex != baseIndex) {
            //Propagate associated joint to the top, otherwise the partial derivative of the last joint is not included in the first joint
            ancestor_ptr = m_traversal.getParentLinkFromLinkIndex(parentLinkIndex);

            while (ancestor_ptr && (ancestor_ptr->getIndex() != baseIndex)) {
                ancestorLink = static_cast<size_t>(ancestor_ptr->getIndex());
                ancestor_T_associated = m_kinDyn.getRelativeTransform(ancestor_ptr->getIndex(), associatedLinkIndex);
                m_childrenForceDerivatives[ancestorLink][associatedJoint] = m_childrenForceDerivatives[ancestorLink][associatedJoint] + (ancestor_T_associated * m_childrenForceDerivatives[associatedLink][associatedJoint]);
                ancestor_ptr = m_traversal.getParentLinkFromLinkIndex(ancestor_ptr->getIndex());
            }

        }
    }

    iDynTree::LinkIndex jLink, parentOfJ;
    iDynTree::SpatialMotionVector sJ;
    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {

        JointInfos& jInfo = m_jointsInfos[j];
        jLink = jInfo.childIndex;
        parentOfJ = jInfo.parentIndex;
        sJ = jInfo.jointPtr->getMotionSubspaceVector(0, jLink, parentOfJ);

        for (size_t z = 0; z < m_jointsInfos.size(); ++z) {
            staticTorquesDerivatives(static_cast<unsigned int>(j), static_cast<unsigned int>(z)) = sJ.dot(m_childrenForceDerivatives[static_cast<size_t>(jLink)][z]);
        }
    }

    return true;
}

bool SharedKinDynComputation::getFreeFloatingMassMatrix(const RobotState &currentState, iDynTree::MatrixDynSize &freeFloatingMassMatrix, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!m_kinDyn.isValid())
        return false;

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (!updateRobotState(currentState))
        return false;

    return m_kinDyn.getFreeFloatingMassMatrix(freeFloatingMassMatrix);
}



