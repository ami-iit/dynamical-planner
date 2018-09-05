/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
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
    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {
        iDynTree::JointIndex jointIndex = static_cast<iDynTree::JointIndex>(j);
        assert(model.isValidJointIndex(jointIndex));
        m_jointsInfos[j].jointPtr = model.getJoint(jointIndex);
        assert(m_jointsInfos[j].jointPtr->getNrOfDOFs() == 1);
        m_jointsInfos[j].childIndex =  m_traversal.getChildLinkIndexFromJointIndex(model, jointIndex);
        m_jointsInfos[j].parentIndex =  m_traversal.getParentLinkIndexFromJointIndex(model, jointIndex);
    }
}

void SharedKinDynComputation::updateChildInformations()
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

    updateChildInformations();

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



