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
#include <mutex>

namespace DynamicalPlanner {
    namespace Private {
        class SharedKinDynComputation;

        typedef struct {
            iDynTree::Transform world_T_base;
            iDynTree::VectorDynSize s;
            iDynTree::Twist base_velocity;
            iDynTree::VectorDynSize s_dot;
        } RobotState;
    }
}

class DynamicalPlanner::Private::SharedKinDynComputation {

    iDynTree::KinDynComputations m_kinDyn;
    std::mutex m_mutex;
    RobotState m_state;
    iDynTree::Vector3 m_gravity;
    bool m_updateNecessary;
    double m_tol;

    bool sameState(const RobotState& other);

public:

    SharedKinDynComputation();

    //Setup

    bool loadRobotModel(const iDynTree::Model& model);

    const iDynTree::Model& model() const;

    bool isValid() const;

    void setGravity(const iDynTree::Vector3& gravity);

    bool setToleranceForUpdate(double tol);

    bool setFloatingBase(const std::string & floatingBaseName);

    //

    const RobotState &currentState() const;

    iDynTree::Position getCenterOfMassPosition(const RobotState &currentState);

    bool getCenterOfMassJacobian(const RobotState &currentState,
                                 iDynTree::MatrixDynSize & comJacobian,
                                 iDynTree::FrameVelocityRepresentation trivialization =
                                    iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

};



#endif // DPLANNER_SHAREDKINDYNCOMPUTATIONS_H
