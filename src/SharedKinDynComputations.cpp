/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/SharedKinDynComputations.h>
#include <private/CheckEqualVector.h>
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

    return m_kinDyn.setFloatingBase(floatingBaseName);
}

const RobotState &SharedKinDynComputation::currentState() const
{
    return m_state;
}

iDynTree::Position SharedKinDynComputation::getCenterOfMassPosition(const RobotState &currentState)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (m_updateNecessary || !sameState(currentState)) {
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        assert(ok);
        m_updateNecessary = false;
        m_state = currentState;
    }

    return m_kinDyn.getCenterOfMassPosition();
}

bool SharedKinDynComputation::getCenterOfMassJacobian(const RobotState &currentState, iDynTree::MatrixDynSize &comJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (m_updateNecessary || !sameState(currentState)) {
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        if (!ok) {
            return false;
        }
        m_updateNecessary = false;
        m_state = currentState;
    }

    return m_kinDyn.getCenterOfMassJacobian(comJacobian);

}

iDynTree::Transform SharedKinDynComputation::getWorldTransform(const RobotState &currentState, std::string frameName)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (m_updateNecessary || !sameState(currentState)) {
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        assert(ok);
        m_updateNecessary = false;
        m_state = currentState;
    }

    return m_kinDyn.getWorldTransform(frameName);
}

iDynTree::Transform SharedKinDynComputation::getWorldTransform(const RobotState &currentState, const iDynTree::FrameIndex frameIndex)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    if (m_updateNecessary || !sameState(currentState)) {
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        assert(ok);
        m_updateNecessary = false;
        m_state = currentState;
    }

    return m_kinDyn.getWorldTransform(frameIndex);
}

bool SharedKinDynComputation::getFrameFreeFloatingJacobian(const RobotState &currentState, const std::string &frameName, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (m_updateNecessary || !sameState(currentState)) {
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        if (!ok) {
            return false;
        }
        m_updateNecessary = false;
        m_state = currentState;
    }

    return m_kinDyn.getFrameFreeFloatingJacobian(frameName, outJacobian);
}

bool SharedKinDynComputation::getFrameFreeFloatingJacobian(const RobotState &currentState, const iDynTree::FrameIndex frameIndex, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_mutex);

    m_kinDyn.setFrameVelocityRepresentation(trivialization);

    if (m_updateNecessary || !sameState(currentState)) {
        bool ok = m_kinDyn.setRobotState(currentState.world_T_base, currentState.s, currentState.base_velocity, currentState.s_dot, m_gravity);
        if (!ok) {
            return false;
        }
        m_updateNecessary = false;
        m_state = currentState;
    }

    return m_kinDyn.getFrameFreeFloatingJacobian(frameIndex, outJacobian);
}


