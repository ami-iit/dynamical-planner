/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SHAREDKINDYNCOMPUTATIONS_H
#define DPLANNER_SHAREDKINDYNCOMPUTATIONS_H

#include <iDynTree/KinDynComputations.h>
#include <iDynTree/Model/Model.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/Twist.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/Model/Traversal.h>
#include <iDynTree/Model/LinkState.h>
#include <iDynTree/Model/FreeFloatingState.h>
#include <mutex>
#include <vector>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class SharedKinDynComputations;
        typedef std::shared_ptr<SharedKinDynComputations> SharedKinDynComputationsPointer;

        typedef struct {
            iDynTree::Transform world_T_base;
            iDynTree::VectorDynSize s;
            iDynTree::Twist base_velocity;
            iDynTree::VectorDynSize s_dot;
        } RobotState;

        typedef struct {
            iDynTree::IJointConstPtr jointPtr;
            iDynTree::SpatialMomentum successorsMomentum, velocityDerivative;
            iDynTree::SpatialMotionVector motionVectorTimesChildVelocity, motionVectorTimesChildAcceleration;
            iDynTree::LinkIndex childIndex, parentIndex;
            iDynTree::Twist childVelocity;
            iDynTree::Transform baseTC;
            iDynTree::SpatialForceVector childStaticForceDerivative;
        } JointInfos;
    }
}

class DynamicalPlanner::Private::SharedKinDynComputations {

    iDynTree::KinDynComputations m_kinDyn;
    std::mutex m_mutex;
    RobotState m_state;
    iDynTree::Vector3 m_gravity;
    std::vector<JointInfos> m_jointsInfos;
    iDynTree::LinkWrenches m_linkStaticWrenches;
    iDynTree::Traversal m_traversal;
    iDynTree::FreeFloatingAcc m_invDynGeneralizedProperAccs;
    iDynTree::Vector3 m_gravityAccInBaseLinkFrame;
    iDynTree::FreeFloatingPos m_pos;
    iDynTree::FreeFloatingVel m_invDynZeroVel;
    iDynTree::LinkVelArray m_invDynZeroLinkVel;
    iDynTree::LinkProperAccArray m_invDynLinkProperAccs;
    iDynTree::FreeFloatingGeneralizedTorques m_generalizedStaticTorques;
    std::vector<std::vector<iDynTree::SpatialForceVector>> m_childrenForceDerivatives;
    std::vector<iDynTree::SpatialForceVector> m_zeroDerivatives;

    bool m_updateNecessary;
    double m_tol;

    bool sameState(const RobotState& other);

    bool updateRobotState(const RobotState &currentState);

    void fillJointsInfo();

    void updateChildBuffersForMomentumDerivative();

    void computeChildStaticForceDerivative(const iDynTree::LinkWrenches &linkStaticForces);

    bool computeStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces);


public:

    SharedKinDynComputations();

    SharedKinDynComputations(const SharedKinDynComputations& other);

    SharedKinDynComputations(const SharedKinDynComputations&& other) = delete;

    //Setup

    bool loadRobotModel(const iDynTree::Model& model);

    const iDynTree::Model& model() const;

    bool isValid() const;

    void setGravity(const iDynTree::Vector3& gravity);

    const iDynTree::Vector3 &gravity() const;

    bool setToleranceForUpdate(double tol);

    double getUpdateTolerance() const;

    bool setFloatingBase(const std::string & floatingBaseName);

    std::string getFloatingBase() const;

    //

    const iDynTree::Traversal &traversal() const;

    const RobotState &currentState() const;

    iDynTree::Position getCenterOfMassPosition(const RobotState &currentState);

    bool getCenterOfMassJacobian(const RobotState &currentState,
                                 iDynTree::MatrixDynSize & comJacobian,
                                 iDynTree::FrameVelocityRepresentation trivialization =
                                    iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    iDynTree::Transform getWorldTransform(const RobotState &currentState, std::string frameName);

    iDynTree::Transform getWorldTransform(const RobotState &currentState, const iDynTree::FrameIndex frameIndex);

    bool getFrameFreeFloatingJacobian(const RobotState &currentState,
                                      const std::string & frameName,
                                      iDynTree::MatrixDynSize & outJacobian,
                                      iDynTree::FrameVelocityRepresentation trivialization =
                                         iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    bool getFrameFreeFloatingJacobian(const RobotState &currentState,
                                      const iDynTree::FrameIndex frameIndex,
                                      iDynTree::MatrixDynSize & outJacobian,
                                      iDynTree::FrameVelocityRepresentation trivialization =
                                         iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    iDynTree::Transform getRelativeTransform(const RobotState &currentState,
                                             const iDynTree::FrameIndex refFrameIndex,
                                             const iDynTree::FrameIndex frameIndex);

    iDynTree::Transform getRelativeTransform(const RobotState &currentState,
                                             const std::string & refFrameName,
                                             const std::string & frameName);

    bool getRelativeJacobian(const RobotState &currentState,
                             const iDynTree::FrameIndex refFrameIndex,
                             const iDynTree::FrameIndex frameIndex,
                             iDynTree::MatrixDynSize & outJacobian,
                             iDynTree::FrameVelocityRepresentation trivialization =
                                iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    iDynTree::Twist getFrameVel(const RobotState &currentState, const std::string & frameName,
                                iDynTree::FrameVelocityRepresentation trivialization =
                                   iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    iDynTree::Twist getFrameVel(const RobotState &currentState, const iDynTree::FrameIndex frameIdx,
                                iDynTree::FrameVelocityRepresentation trivialization =
                                   iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    bool getFrameVelJointsDerivative(const RobotState &currentState, const iDynTree::FrameIndex frameIdx,
                                     iDynTree::MatrixDynSize& velocityDerivative); //Implemented only for BODY_REPRESENTATION

    iDynTree::SpatialMomentum getLinearAngularMomentum(const RobotState &currentState,
                                                       iDynTree::FrameVelocityRepresentation trivialization =
                                                           iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    bool getLinearAngularMomentumJacobian(const RobotState &currentState, iDynTree::MatrixDynSize &linAngMomentumJacobian,
                                          iDynTree::FrameVelocityRepresentation trivialization =
                                              iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

    bool getLinearAngularMomentumJointsDerivative(const RobotState &currentState, iDynTree::MatrixDynSize &linAngMomentumDerivative); //Implemented only for BODY_REPRESENTATION

    bool getStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces,
                         iDynTree::FreeFloatingGeneralizedTorques &generalizedStaticTorques, iDynTree::LinkWrenches &linkStaticForces); //The external forces are expected in an BODY_REPRESENTATION. The base generalized torque and the linkStaticForces are expressed in BODY_REPRESENTATION

    bool getStaticForces(const RobotState &currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces,
                         iDynTree::FreeFloatingGeneralizedTorques &generalizedStaticTorques); //The external forces are expected in an BODY_REPRESENTATION. The base generalized torque and the linkStaticForces are expressed in BODY_REPRESENTATION

    bool getStaticForcesJointsDerivative(const RobotState& currentState, const iDynTree::LinkNetExternalWrenches &linkExtForces,
                                         iDynTree::MatrixDynSize &staticTorquesDerivatives); //The external forces are expected in an BODY_REPRESENTATION. The base generalized torque and the linkStaticForces are expressed in BODY_REPRESENTATION

    bool getFreeFloatingMassMatrix(const RobotState& currentState, iDynTree::MatrixDynSize & freeFloatingMassMatrix,
                                   iDynTree::FrameVelocityRepresentation trivialization =
            iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);



};



#endif // DPLANNER_SHAREDKINDYNCOMPUTATIONS_H
