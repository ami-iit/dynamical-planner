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

void SharedKinDynComputation::getChildSpatialMotionVectors()
{
    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {

        m_jointsInfos[j].childMotionVector = m_jointsInfos[j].jointPtr->getMotionSubspaceVector(0,
                                                                                                m_jointsInfos[j].childIndex,
                                                                                                m_jointsInfos[j].parentIndex).cross(m_kinDyn.getFrameVel(m_jointsInfos[j].childIndex));
    }
}

void SharedKinDynComputation::resetVisits()
{
    for (unsigned int j = 0; j < m_jointsInfos.size(); ++j) {
        m_jointsInfos[j].alreadyVisited = false;
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

bool SharedKinDynComputation::getLinearAngularMomentumJointsDerivative(const RobotState &currentState, iDynTree::MatrixDynSize &linAngMomentumDerivative)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (!m_kinDyn.isValid())
        return false;

    m_kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    if (!updateRobotState(currentState))
        return false;

    resetVisits();
    getChildSpatialMotionVectors();

    const iDynTree::Model& model = m_kinDyn.model();
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    iDynTree::IJointConstPtr jointPtr;
    iDynTree::LinkConstPtr linkPtr;
    size_t jointIndex;
    for (size_t l = 0; l < model.getNrOfLinks(); ++l) {

        iDynTree::LinkIndex linkIndex = static_cast<iDynTree::LinkIndex>(l);
        assert(model.isValidLinkIndex(linkIndex));
        linkPtr = model.getLink(linkIndex);

        while (linkIndex != baseIndex) {
            jointPtr = m_traversal.getParentJointFromLinkIndex(linkIndex);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());
            if (!m_jointsInfos[jointIndex].alreadyVisited) {
                m_jointsInfos[jointIndex].successorsMomentum = m_kinDyn.getRelativeTransform(m_jointsInfos[jointIndex].childIndex, linkIndex) *
                        linkPtr->getInertia() * m_kinDyn.getFrameVel(linkIndex);

                m_jointsInfos[jointIndex].velocityDerivative = m_kinDyn.getRelativeTransform(baseIndex, linkIndex) * linkPtr->getInertia() *
                        (m_kinDyn.getRelativeTransform(linkIndex, m_jointsInfos[jointIndex].childIndex) * m_jointsInfos[jointIndex].childMotionVector);


                m_jointsInfos[jointIndex].alreadyVisited = true;
            } else {
                m_jointsInfos[jointIndex].successorsMomentum = m_jointsInfos[jointIndex].successorsMomentum +
                        m_kinDyn.getRelativeTransform(m_jointsInfos[jointIndex].childIndex, linkIndex) *
                        linkPtr->getInertia() * m_kinDyn.getFrameVel(linkIndex);

                m_jointsInfos[jointIndex].velocityDerivative = m_jointsInfos[jointIndex].velocityDerivative + m_kinDyn.getRelativeTransform(baseIndex, linkIndex) * linkPtr->getInertia() *
                        (m_kinDyn.getRelativeTransform(linkIndex, m_jointsInfos[jointIndex].childIndex) * m_jointsInfos[jointIndex].childMotionVector);

            }

            linkIndex = m_traversal.getParentLinkFromLinkIndex(linkIndex)->getIndex();

        }

    }

    linAngMomentumDerivative.resize(6, static_cast<unsigned int>(m_jointsInfos.size()));
    iDynTree::iDynTreeEigenMatrixMap derivativeMap = iDynTree::toEigen(linAngMomentumDerivative);
    iDynTree::SpatialMomentum colBuffer;

    for (size_t j = 0; j < m_jointsInfos.size(); ++j) {
        assert(m_jointsInfos[j].alreadyVisited);

        colBuffer = m_kinDyn.getRelativeTransform(baseIndex, m_jointsInfos[j].childIndex) *
                m_jointsInfos[j].jointPtr->getMotionSubspaceVector(0,
                                                                   m_jointsInfos[j].childIndex,
                                                                   m_jointsInfos[j].parentIndex).cross(m_jointsInfos[j].successorsMomentum) +
                m_jointsInfos[j].velocityDerivative;

        derivativeMap.col(static_cast<Eigen::Index>(j)) = iDynTree::toEigen(colBuffer);
    }

    return true;
}



