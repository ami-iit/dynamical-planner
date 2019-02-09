/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>
#include <iDynTree/Model/ForwardKinematics.h>
#include <iDynTree/Model/Dynamics.h>
#include <DynamicalPlannerPrivate/CheckEqualVector.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>

using namespace DynamicalPlanner::Private;

class SharedKinDynComputations::Data {
public:
    iDynTree::KinDynComputations kinDyn;
    std::mutex mutex;
    RobotState state;
    iDynTree::Vector3 gravity;
    std::vector<JointInfos> jointsInfos;
    iDynTree::LinkWrenches linkStaticWrenches;
    iDynTree::Traversal traversal;
    iDynTree::FreeFloatingAcc invDynGeneralizedProperAccs;
    iDynTree::Vector3 gravityAccInBaseLinkFrame;
    iDynTree::FreeFloatingPos pos;
    iDynTree::FreeFloatingVel invDynZeroVel;
    iDynTree::LinkVelArray invDynZeroLinkVel;
    iDynTree::LinkProperAccArray invDynLinkProperAccs;
    iDynTree::FreeFloatingGeneralizedTorques generalizedStaticTorques;
    std::vector<std::vector<iDynTree::SpatialForceVector>> childrenForceDerivatives;
    std::vector<iDynTree::SpatialForceVector> zeroDerivatives;
    iDynTree::Transform baseTransform;
    iDynTree::Rotation baseRotation_iDyn;
    iDynTree::Position basePosition;

    levi::Variable quaternion = levi::Variable(4, "q");
    levi::Expression quaternionNormalized;
    levi::Expression baseRotation;

    bool updateNecessary;
    double tol;
};

bool SharedKinDynComputations::sameState(const RobotState &other)
{
    if (VectorsAreEqual(other.base_position, m_data->state.base_position, m_data->tol)
            && VectorsAreEqual(other.base_quaternion, m_data->state.base_quaternion, m_data->tol)
            && VectorsAreEqual(other.s, m_data->state.s, m_data->tol)
            && VectorsAreEqual(other.base_velocity, m_data->state.base_velocity, m_data->tol)
            && VectorsAreEqual(other.s_dot, m_data->state.s_dot, m_data->tol)) {
        return true;
    }
    return false;
}

bool SharedKinDynComputations::updateRobotState(const RobotState &currentState)
{
    if (m_data->updateNecessary || !sameState(currentState)) {

        m_data->quaternion = iDynTree::toEigen(currentState.base_quaternion);
        iDynTree::toEigen(m_data->baseRotation_iDyn) = m_data->baseRotation.evaluate();
        iDynTree::toEigen(m_data->basePosition) = iDynTree::toEigen(currentState.base_position);
        m_data->baseTransform.setRotation(m_data->baseRotation_iDyn);
        m_data->baseTransform.setPosition(m_data->basePosition);

        m_data->kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION); //The base_velocity saved in the robot state is supposed to be in body frame
        bool ok = m_data->kinDyn.setRobotState(m_data->baseTransform, currentState.s, currentState.base_velocity, currentState.s_dot, m_data->gravity);
        if (!ok) {
            return false;
        }
        m_data->updateNecessary = false;
        m_data->state = currentState;
    }
    return true;
}

void SharedKinDynComputations::fillJointsInfo()
{
    const iDynTree::Model& model = m_data->kinDyn.model();
    m_data->childrenForceDerivatives.resize(model.getNrOfLinks());
    for (size_t j = 0; j < m_data->jointsInfos.size(); ++j) {
        iDynTree::JointIndex jointIndex = static_cast<iDynTree::JointIndex>(j);
        assert(model.isValidJointIndex(jointIndex));
        m_data->jointsInfos[j].jointPtr = model.getJoint(jointIndex);
        assert(m_data->jointsInfos[j].jointPtr->getNrOfDOFs() == 1);
        m_data->jointsInfos[j].childIndex =  m_data->traversal.getChildLinkIndexFromJointIndex(model, jointIndex);
        m_data->jointsInfos[j].parentIndex =  m_data->traversal.getParentLinkIndexFromJointIndex(model, jointIndex);
        size_t childIndex = static_cast<size_t>(m_data->jointsInfos[j].childIndex);
        m_data->childrenForceDerivatives[childIndex].resize(m_data->jointsInfos.size(), iDynTree::SpatialForceVector::Zero());
    }
}

void SharedKinDynComputations::updateChildBuffersForMomentumDerivative()
{
    iDynTree::LinkIndex baseIndex = m_data->kinDyn.model().getLinkIndex(m_data->kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    for (size_t j = 0; j < m_data->jointsInfos.size(); ++j) {

        m_data->jointsInfos[j].motionVectorTimesChildVelocity = m_data->jointsInfos[j].jointPtr->getMotionSubspaceVector(0,
                                                                                                m_data->jointsInfos[j].childIndex,
                                                                                                m_data->jointsInfos[j].parentIndex).cross(
                    -m_data->kinDyn.getFrameVel(m_data->jointsInfos[j].childIndex));

        m_data->jointsInfos[j].childVelocity = m_data->kinDyn.getFrameVel(m_data->jointsInfos[j].childIndex);

        m_data->jointsInfos[j].baseTC = m_data->kinDyn.getRelativeTransform(baseIndex, m_data->jointsInfos[j].childIndex);
    }
}

void SharedKinDynComputations::computeChildStaticForceDerivative(const iDynTree::LinkWrenches &linkStaticForces)
{

    iDynTree::LinkIndex childLink, parentLink;
    size_t childIndex;
    for (size_t j = 0; j < m_data->jointsInfos.size(); ++j) {

        childLink = m_data->jointsInfos[j].childIndex;
        parentLink = m_data->jointsInfos[j].parentIndex;

        JointInfos& jInfo = m_data->jointsInfos[j];

        jInfo.childStaticForceDerivative = jInfo.jointPtr->getMotionSubspaceVector(0,
                                                                                   childLink,
                                                                                   parentLink).cross(
                    linkStaticForces(childLink));

        jInfo.motionVectorTimesChildAcceleration = jInfo.jointPtr->getMotionSubspaceVector(0,
                                                                                           childLink,
                                                                                           parentLink).cross(
                    m_data->invDynLinkProperAccs(childLink));

        childIndex = static_cast<size_t>(childLink);
        m_data->childrenForceDerivatives[childIndex] = m_data->zeroDerivatives;
    }
}

bool SharedKinDynComputations::computeStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces)
{
    iDynTree::Transform baseTransform;
    iDynTree::Vector4 inputQuaternion;
    iDynTree::Rotation baseRotation;
    iDynTree::Position basePosition;

    iDynTree::toEigen(inputQuaternion) = iDynTree::toEigen(currentState.base_quaternion).normalized();
    baseRotation.fromQuaternion(inputQuaternion);
    iDynTree::toEigen(basePosition) = iDynTree::toEigen(currentState.base_position);
    baseTransform.setPosition(basePosition);
    baseTransform.setRotation(baseRotation);

    iDynTree::toEigen(m_data->gravityAccInBaseLinkFrame) =
            iDynTree::toEigen(baseRotation.inverse())*iDynTree::toEigen(m_data->gravity);

    // Clear input buffers that need to be cleared
    m_data->invDynGeneralizedProperAccs.baseAcc().zero();
    iDynTree::toEigen(m_data->invDynGeneralizedProperAccs.baseAcc().getLinearVec3()) = - iDynTree::toEigen(m_data->gravityAccInBaseLinkFrame);
    m_data->invDynGeneralizedProperAccs.jointAcc().zero();

    m_data->pos.worldBasePos() = baseTransform;
    iDynTree::toEigen(m_data->pos.jointPos()) = iDynTree::toEigen(currentState.s);

    // Run inverse dynamics
    bool ok = iDynTree::ForwardAccKinematics(m_data->kinDyn.model(),
                                             m_data->traversal,
                                             m_data->pos,
                                             m_data->invDynZeroVel,
                                             m_data->invDynGeneralizedProperAccs,
                                             m_data->invDynZeroLinkVel,
                                             m_data->invDynLinkProperAccs);

    ok = ok && iDynTree::RNEADynamicPhase(m_data->kinDyn.model(),
                                          m_data->traversal,
                                          m_data->pos.jointPos(),
                                          m_data->invDynZeroLinkVel,
                                          m_data->invDynLinkProperAccs,
                                          linkExtForces,
                                          m_data->linkStaticWrenches,
                                          m_data->generalizedStaticTorques);
    return ok;
}

SharedKinDynComputations::SharedKinDynComputations()
    : m_data(std::make_unique<Data>())
{
    m_data->state.base_quaternion.zero();
    m_data->state.base_quaternion(0) = 1.0;
    m_data->state.base_position.zero();
    m_data->state.base_velocity.zero();
    m_data->gravity.zero();
    m_data->gravity(2) = -9.81;
    m_data->updateNecessary = true;
    m_data->tol = iDynTree::DEFAULT_TOL;

    m_data->quaternionNormalized = m_data->quaternion/(m_data->quaternion.transpose() * m_data->quaternion).pow(0.5);
    levi::Expression skewQuaternion = m_data->quaternionNormalized.block(1,0,3,1).skew();
    levi::Expression twoSkewQuaternion = 2.0 * skewQuaternion;
    m_data->baseRotation = levi::Identity(3,3) + m_data->quaternionNormalized(0,0) * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;

}

SharedKinDynComputations::SharedKinDynComputations(const DynamicalPlanner::Private::SharedKinDynComputations &other)
    : m_data(std::make_unique<Data>())
{
    assert(other.isValid());

    m_data->state.base_quaternion.zero();
    m_data->state.base_quaternion(0) = 1.0;
    m_data->state.base_position.zero();
    m_data->state.base_velocity.zero();
    m_data->gravity = other.gravity();
    m_data->updateNecessary = true;

    loadRobotModel(other.model());

    setToleranceForUpdate(other.getUpdateTolerance());

    setFloatingBase(other.getFloatingBase());

    assert(isValid());

    m_data->quaternionNormalized = m_data->quaternion/(m_data->quaternion.transpose() * m_data->quaternion).pow(0.5);
    levi::Expression skewQuaternion = m_data->quaternionNormalized.block(1,0,3,1).skew();
    levi::Expression twoSkewQuaternion = 2.0 * skewQuaternion;
    m_data->baseRotation = levi::Identity(3,3) + m_data->quaternionNormalized(0,0) * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;
}

SharedKinDynComputations::~SharedKinDynComputations()
{ }

bool SharedKinDynComputations::loadRobotModel(const iDynTree::Model &model)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);
    bool ok = m_data->kinDyn.loadRobotModel(model);

    if (ok) {
        m_data->state.s.resize(m_data->kinDyn.getNrOfDegreesOfFreedom());
        m_data->state.s.zero();
        m_data->state.s_dot.resize(m_data->kinDyn.getNrOfDegreesOfFreedom());
        m_data->state.s_dot.zero();
        m_data->jointsInfos.resize(m_data->kinDyn.getNrOfDegreesOfFreedom());
        bool okTraversal = m_data->kinDyn.model().computeFullTreeTraversal(m_data->traversal);
        if (!okTraversal)
            return false;
        fillJointsInfo();
    }

    m_data->linkStaticWrenches.resize(model);
    m_data->invDynGeneralizedProperAccs.resize(model);
    m_data->pos.resize(model);
    m_data->invDynZeroVel.resize(model);
    m_data->invDynZeroVel.baseVel().zero();
    m_data->invDynZeroVel.jointVel().zero();
    m_data->invDynZeroLinkVel.resize(model);
    m_data->invDynLinkProperAccs.resize(model);

    for(iDynTree::LinkIndex lnkIdx = 0; lnkIdx < static_cast<iDynTree::LinkIndex>(model.getNrOfLinks()); lnkIdx++)
    {
        m_data->invDynZeroLinkVel(lnkIdx).zero();
    }

    m_data->generalizedStaticTorques.resize(model);
    m_data->zeroDerivatives.resize(m_data->kinDyn.getNrOfDegreesOfFreedom(), iDynTree::SpatialForceVector::Zero());

    return ok;
}

const iDynTree::Model &SharedKinDynComputations::model() const
{
    return m_data->kinDyn.model();
}

bool SharedKinDynComputations::isValid() const
{
    return m_data->kinDyn.isValid();
}

void SharedKinDynComputations::setGravity(const iDynTree::Vector3 &gravity)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);
    if (!VectorsAreEqual(gravity, m_data->gravity, m_data->tol))
        m_data->updateNecessary = true;

    m_data->gravity = gravity;
}

const iDynTree::Vector3 &SharedKinDynComputations::gravity() const
{
    return m_data->gravity;
}

bool SharedKinDynComputations::setToleranceForUpdate(double tol)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (tol < 0) {
        return false;
    }

    m_data->tol = tol;

    return true;
}

double SharedKinDynComputations::getUpdateTolerance() const
{
    return m_data->tol;
}

bool SharedKinDynComputations::setFloatingBase(const std::string &floatingBaseName)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    m_data->updateNecessary = true;

    if (!m_data->kinDyn.setFloatingBase(floatingBaseName)) {
        return false;
    }

    bool okTraversal = m_data->kinDyn.model().computeFullTreeTraversal(m_data->traversal, m_data->kinDyn.model().getLinkIndex(m_data->kinDyn.getFloatingBase()));
    if (!okTraversal)
        return false;

    fillJointsInfo();

    return true;
}

std::string SharedKinDynComputations::getFloatingBase() const
{
    return m_data->kinDyn.getFloatingBase();
}

const iDynTree::Traversal &SharedKinDynComputations::traversal() const
{
    return m_data->traversal;
}

const RobotState &SharedKinDynComputations::currentState() const
{
    return m_data->state;
}

iDynTree::Position SharedKinDynComputations::getCenterOfMassPosition(const RobotState &currentState)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->kinDyn.getCenterOfMassPosition();
}

bool SharedKinDynComputations::getCenterOfMassJacobian(const RobotState &currentState, iDynTree::MatrixDynSize &comJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getCenterOfMassJacobian(comJacobian);

}

iDynTree::Transform SharedKinDynComputations::getWorldTransform(const RobotState &currentState, std::string frameName)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->kinDyn.getWorldTransform(frameName);
}

iDynTree::Transform SharedKinDynComputations::getWorldTransform(const RobotState &currentState, const iDynTree::FrameIndex frameIndex)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->kinDyn.getWorldTransform(frameIndex);
}

const iDynTree::Transform &SharedKinDynComputations::getBaseTransform(const RobotState &currentState)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->baseTransform;
}

levi::Expression& SharedKinDynComputations::baseRotation(const RobotState &currentState)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->baseRotation;
}

levi::Variable &SharedKinDynComputations::baseQuaternion(const RobotState &currentState)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->quaternion;
}

bool SharedKinDynComputations::getFrameFreeFloatingJacobian(const RobotState &currentState, const std::string &frameName, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getFrameFreeFloatingJacobian(frameName, outJacobian);
}

bool SharedKinDynComputations::getFrameFreeFloatingJacobian(const RobotState &currentState, const iDynTree::FrameIndex frameIndex, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getFrameFreeFloatingJacobian(frameIndex, outJacobian);
}

iDynTree::Transform SharedKinDynComputations::getRelativeTransform(const RobotState &currentState,
                                                                  const iDynTree::FrameIndex refFrameIndex,
                                                                  const iDynTree::FrameIndex frameIndex)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->kinDyn.getRelativeTransform(refFrameIndex, frameIndex);
}

iDynTree::Transform SharedKinDynComputations::getRelativeTransform(const RobotState &currentState,
                                                                  const std::string &refFrameName,
                                                                  const std::string &frameName)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    return m_data->kinDyn.getRelativeTransform(refFrameName, frameName);
}

bool SharedKinDynComputations::getRelativeJacobian(const RobotState &currentState, const iDynTree::FrameIndex refFrameIndex, const iDynTree::FrameIndex frameIndex, iDynTree::MatrixDynSize &outJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getRelativeJacobian(refFrameIndex, frameIndex, outJacobian);
}

iDynTree::Twist SharedKinDynComputations::getFrameVel(const RobotState &currentState, const std::string &frameName, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getFrameVel(frameName);
}

iDynTree::Twist SharedKinDynComputations::getFrameVel(const RobotState &currentState, const iDynTree::FrameIndex frameIdx, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getFrameVel(frameIdx);
}

bool SharedKinDynComputations::getFrameVelJointsDerivative(const RobotState &currentState, const iDynTree::FrameIndex frameIdx, iDynTree::MatrixDynSize &velocityDerivative)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!m_data->kinDyn.isValid())
        return false;

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    velocityDerivative.resize(6, static_cast<unsigned int>(m_data->jointsInfos.size()));
    velocityDerivative.zero();
    iDynTree::iDynTreeEigenMatrixMap derivativeMap = iDynTree::toEigen(velocityDerivative);

    const iDynTree::Model& model = m_data->kinDyn.model();
    iDynTree::LinkIndex linkIndex = model.getFrameLink(frameIdx);
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_data->kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    if (linkIndex == iDynTree::LINK_INVALID_INDEX) {
        return false;
    }

    iDynTree::IJointConstPtr jointPtr;
    size_t jointIndex;

    iDynTree::LinkIndex visitedLink = linkIndex;
    iDynTree::Twist jointDerivative, childVelocity;
    while (visitedLink != baseIndex) {
        jointPtr = m_data->traversal.getParentJointFromLinkIndex(visitedLink);
        jointIndex = static_cast<size_t>(jointPtr->getIndex());

        childVelocity = m_data->kinDyn.getFrameVel(m_data->jointsInfos[jointIndex].childIndex);
        jointDerivative = m_data->kinDyn.getRelativeTransform(frameIdx, m_data->jointsInfos[jointIndex].childIndex) *
                m_data->jointsInfos[jointIndex].jointPtr->getMotionSubspaceVector(0,
                                                                            m_data->jointsInfos[jointIndex].childIndex,
                                                                            m_data->jointsInfos[jointIndex].parentIndex).cross(-childVelocity);

        derivativeMap.col(static_cast<Eigen::Index>(jointIndex)) = iDynTree::toEigen(jointDerivative);
        visitedLink = m_data->traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
    }

    return true;
}

iDynTree::SpatialMomentum SharedKinDynComputations::getLinearAngularMomentum(const RobotState &currentState, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    bool ok = updateRobotState(currentState);
    assert(ok);

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getLinearAngularMomentum();
}

bool SharedKinDynComputations::getLinearAngularMomentumJacobian(const RobotState &currentState, iDynTree::MatrixDynSize &linAngMomentumJacobian, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getLinearAngularMomentumJacobian(linAngMomentumJacobian);
}

bool SharedKinDynComputations::getLinearAngularMomentumJointsDerivative(const RobotState &currentState,
                                                                       iDynTree::MatrixDynSize &linAngMomentumDerivative)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!m_data->kinDyn.isValid())
        return false;

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    linAngMomentumDerivative.resize(6, static_cast<unsigned int>(m_data->jointsInfos.size()));
    linAngMomentumDerivative.zero();

    updateChildBuffersForMomentumDerivative();

    iDynTree::iDynTreeEigenMatrixMap derivativeMap = iDynTree::toEigen(linAngMomentumDerivative);

    const iDynTree::Model& model = m_data->kinDyn.model();
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_data->kinDyn.getFloatingBase());
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
        linkMomentum = linkPtr->getInertia() * m_data->kinDyn.getFrameVel(linkIndex);
        b_T_link = m_data->kinDyn.getRelativeTransform(baseIndex, linkIndex);

        visitedLink = linkIndex;
        while(visitedLink != baseIndex) {
            jointPtr = m_data->traversal.getParentJointFromLinkIndex(visitedLink);
            jointIndex = static_cast<size_t>(jointPtr->getIndex());

            childLink = m_data->jointsInfos[jointIndex].childIndex;
            parentLink = m_data->jointsInfos[jointIndex].parentIndex;

            childVelocity = m_data->jointsInfos[jointIndex].childVelocity;
            linkMomentumInChildFrame = m_data->kinDyn.getRelativeTransform(childLink, linkIndex) * linkMomentum;

            transformDerivative = m_data->jointsInfos[jointIndex].baseTC *
                    jointPtr->getMotionSubspaceVector(0,
                                                      childLink,
                                                      parentLink).cross(linkMomentumInChildFrame);

            velocityDerivative = b_T_link * (linkPtr->getInertia() * (m_data->kinDyn.getRelativeTransform(linkIndex, childLink) *
                                                                      m_data->jointsInfos[jointIndex].motionVectorTimesChildVelocity));

            jointMomentumDerivative = transformDerivative + velocityDerivative;

            derivativeMap.col(static_cast<Eigen::Index>(jointIndex)) += iDynTree::toEigen(jointMomentumDerivative);

            visitedLink = m_data->traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
        }
    }
    return true;
}

bool SharedKinDynComputations::getStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces, iDynTree::FreeFloatingGeneralizedTorques& generalizedStaticTorques, iDynTree::LinkWrenches& linkStaticForces)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!computeStaticForces(currentState, linkExtForces)) {
        return false;
    }

    linkStaticForces = m_data->linkStaticWrenches;
    generalizedStaticTorques = m_data->generalizedStaticTorques;

    return true;
}

bool SharedKinDynComputations::getStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces, iDynTree::FreeFloatingGeneralizedTorques &generalizedStaticTorques)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!computeStaticForces(currentState, linkExtForces)) {
        return false;
    }

    generalizedStaticTorques = m_data->generalizedStaticTorques;

    return true;
}


bool SharedKinDynComputations::getStaticForcesJointsDerivative(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces, iDynTree::MatrixDynSize &staticTorquesDerivatives)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!m_data->kinDyn.isValid())
        return false;

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(iDynTree::FrameVelocityRepresentation::BODY_FIXED_REPRESENTATION);

    if (!computeStaticForces(currentState, linkExtForces)) {
        return false;
    }

    computeChildStaticForceDerivative(m_data->linkStaticWrenches);

    staticTorquesDerivatives.resize(static_cast<unsigned int>(m_data->jointsInfos.size()), static_cast<unsigned int>(m_data->jointsInfos.size()));
    staticTorquesDerivatives.zero();

    const iDynTree::Model& model = m_data->kinDyn.model();
    iDynTree::LinkIndex baseIndex = model.getLinkIndex(m_data->kinDyn.getFloatingBase());
    assert(baseIndex != iDynTree::LINK_INVALID_INDEX);

    iDynTree::LinkIndex parentLinkIndex, associatedLinkIndex;
    iDynTree::LinkConstPtr ancestor_ptr;
    iDynTree::TraversalIndex elIndex;
    iDynTree::IJointConstPtr visitedJoint, associatedJoint_ptr;
    iDynTree::SpatialInertia linkInertia;
    size_t visitedJointIndex, associatedJoint, parentLink, associatedLink, ancestorLink;
    iDynTree::Transform l_T_c, p_T_c, ancestor_T_associated;

    for (unsigned int el = m_data->traversal.getNrOfVisitedLinks() -1; el > 0; el--) {

        elIndex = static_cast<iDynTree::TraversalIndex>(el);
        associatedJoint_ptr = m_data->traversal.getParentJoint(elIndex);
        associatedJoint = static_cast<size_t>(associatedJoint_ptr->getIndex());
        associatedLinkIndex = m_data->traversal.getLink(elIndex)->getIndex();
        linkInertia = model.getLink(associatedLinkIndex)->getInertia();
        associatedLink = static_cast<size_t>(associatedLinkIndex);

        parentLinkIndex = m_data->traversal.getParentLink(elIndex)->getIndex();
        parentLink = static_cast<size_t>(parentLinkIndex);

        p_T_c = m_data->kinDyn.getRelativeTransform(parentLinkIndex, associatedLinkIndex);

        m_data->childrenForceDerivatives[associatedLink][associatedJoint] = m_data->childrenForceDerivatives[associatedLink][associatedJoint] + m_data->jointsInfos[associatedJoint].childStaticForceDerivative;

        // Propagate joints which are before the joint
        visitedJoint = associatedJoint_ptr;
        while (visitedJoint) {
            visitedJointIndex = static_cast<size_t>(visitedJoint->getIndex());

            l_T_c = m_data->kinDyn.getRelativeTransform(associatedLinkIndex, m_data->jointsInfos[visitedJointIndex].childIndex);
            m_data->childrenForceDerivatives[associatedLink][visitedJointIndex] = m_data->childrenForceDerivatives[associatedLink][visitedJointIndex] - linkInertia * (l_T_c * m_data->jointsInfos[visitedJointIndex].motionVectorTimesChildAcceleration);

            if (parentLinkIndex != baseIndex) {
                m_data->childrenForceDerivatives[parentLink][visitedJointIndex] = m_data->childrenForceDerivatives[parentLink][visitedJointIndex] + (p_T_c * m_data->childrenForceDerivatives[associatedLink][visitedJointIndex]);
            }
            visitedJoint = m_data->traversal.getParentJointFromLinkIndex(m_data->jointsInfos[visitedJointIndex].parentIndex);
        }

        if (parentLinkIndex != baseIndex) {
            //Propagate associated joint to the top, otherwise the partial derivative of the last joint is not included in the first joint
            ancestor_ptr = m_data->traversal.getParentLinkFromLinkIndex(parentLinkIndex);

            while (ancestor_ptr && (ancestor_ptr->getIndex() != baseIndex)) {
                ancestorLink = static_cast<size_t>(ancestor_ptr->getIndex());
                ancestor_T_associated = m_data->kinDyn.getRelativeTransform(ancestor_ptr->getIndex(), associatedLinkIndex);
                m_data->childrenForceDerivatives[ancestorLink][associatedJoint] = m_data->childrenForceDerivatives[ancestorLink][associatedJoint] + (ancestor_T_associated * m_data->childrenForceDerivatives[associatedLink][associatedJoint]);
                ancestor_ptr = m_data->traversal.getParentLinkFromLinkIndex(ancestor_ptr->getIndex());
            }

        }
    }

    iDynTree::LinkIndex jLink, parentOfJ;
    iDynTree::SpatialMotionVector sJ;
    for (size_t j = 0; j < m_data->jointsInfos.size(); ++j) {

        JointInfos& jInfo = m_data->jointsInfos[j];
        jLink = jInfo.childIndex;
        parentOfJ = jInfo.parentIndex;
        sJ = jInfo.jointPtr->getMotionSubspaceVector(0, jLink, parentOfJ);

        for (size_t z = 0; z < m_data->jointsInfos.size(); ++z) {
            staticTorquesDerivatives(static_cast<unsigned int>(j), static_cast<unsigned int>(z)) = sJ.dot(m_data->childrenForceDerivatives[static_cast<size_t>(jLink)][z]);
        }
    }

    return true;
}

bool SharedKinDynComputations::getFreeFloatingMassMatrix(const RobotState &currentState, iDynTree::MatrixDynSize &freeFloatingMassMatrix, iDynTree::FrameVelocityRepresentation trivialization)
{
    std::lock_guard<std::mutex> guard(m_data->mutex);

    if (!m_data->kinDyn.isValid())
        return false;

    if (!updateRobotState(currentState))
        return false;

    m_data->kinDyn.setFrameVelocityRepresentation(trivialization);

    return m_data->kinDyn.getFreeFloatingMassMatrix(freeFloatingMassMatrix);
}



