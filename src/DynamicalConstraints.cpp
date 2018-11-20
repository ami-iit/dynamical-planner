/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/DynamicalConstraints.h>
#include <DynamicalPlannerPrivate/QuaternionUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/Utils.h>
#include <cassert>
#include <iostream>
#include <regex>

using namespace DynamicalPlanner::Private;

class DynamicalConstraints::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;
    VariablesLabeller dynamics;
    double totalMass;
    iDynTree::Vector6 gravityVector;

    iDynTree::Position comPosition;
    iDynTree::MatrixDynSize comjacobianBuffer;


    iDynTree::Rotation baseRotation;
    iDynTree::Position basePosition;
    iDynTree::Vector4 baseQuaternion, baseQuaternionNormalized, baseQuaternionVelocity;
    RobotState robotState;
    std::shared_ptr<SharedKinDynComputations> sharedKinDyn;
    std::shared_ptr<TimelySharedKinDynComputations> timedSharedKinDyn;

    HyperbolicTangent activationXY;

    typedef struct {
        std::vector<iDynTree::IndexRange> positionPoints, forcePoints, velocityControlPoints, forceControlPoints;
    } FootRanges;
    FootRanges leftRanges, rightRanges;
    iDynTree::IndexRange momentumRange, comPositionRange, basePositionRange, baseQuaternionRange, jointsPositionRange, jointsVelocityRange;
//    iDynTree::IndexRange baseVelocityRange;
    iDynTree::IndexRange baseLinearVelocityRange, baseQuaternionDerivativeRange;

    //iDynTree::MatrixDynSize stateJacobianBuffer, controlJacobianBuffer;
    iDynTree::optimalcontrol::SparsityStructure stateSparsity, controlSparsity;

    void checkFootVariables(const std::string& footName, size_t numberOfPoints, FootRanges& foot) {
        foot.positionPoints.resize(numberOfPoints);
        foot.forcePoints.resize(numberOfPoints);
        foot.velocityControlPoints.resize(numberOfPoints);
        foot.forceControlPoints.resize(numberOfPoints);

        iDynTree::IndexRange obtainedRange;
        for (size_t i= 0; i < numberOfPoints; ++i) {
            obtainedRange = stateVariables.getIndexRange(footName + "ForcePoint" + std::to_string(i));

            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "ForcePoint" + std::to_string(i) << " not available among state variables." << std::endl;
                assert(false);
            }
            foot.forcePoints[i] = obtainedRange;

            obtainedRange = stateVariables.getIndexRange(footName + "PositionPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "PositionPoint" + std::to_string(i) << " not available among state variables." << std::endl;
                assert(false);
            }
            foot.positionPoints[i] = obtainedRange;

            obtainedRange = controlVariables.getIndexRange(footName + "ForceControlPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "ForceControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
                assert(false);
            }
            foot.forceControlPoints[i] = obtainedRange;

            obtainedRange = controlVariables.getIndexRange(footName + "VelocityControlPoint" + std::to_string(i));
            if (!obtainedRange.isValid()){
                std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << footName + "VelocityControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
                assert(false);
            }
            foot.velocityControlPoints[i] = obtainedRange;
        }
    }


    void computeFootRelatedDynamics(FootRanges& foot) {
        Eigen::Vector3d distance, appliedForce;

        for (size_t i = 0; i < foot.positionPoints.size(); ++i) {
//Span operator = does not copy content!

            iDynTree::toEigen(dynamics(foot.forcePoints[i])) = iDynTree::toEigen(controlVariables(foot.forceControlPoints[i]));
            double deltaXY = activationXY.eval(stateVariables(foot.positionPoints[i])(2));
            iDynTree::toEigen(dynamics(foot.positionPoints[i])).topRows<2>() = deltaXY * iDynTree::toEigen(controlVariables(foot.velocityControlPoints[i])).topRows<2>();
            dynamics(foot.positionPoints[i])(2) = controlVariables(foot.velocityControlPoints[i])(2);

            iDynTree::toEigen(dynamics(momentumRange)).topRows<3>() +=  iDynTree::toEigen(stateVariables(foot.forcePoints[i]));

            distance = iDynTree::toEigen(stateVariables(foot.positionPoints[i])) - iDynTree::toEigen(comPosition);
            appliedForce = iDynTree::toEigen(stateVariables(foot.forcePoints[i]));
            iDynTree::toEigen(dynamics(momentumRange)).bottomRows<3>() += distance.cross(appliedForce);
        }
    }

    void computeFootRelatedStateJacobian(const FootRanges& foot, iDynTree::iDynTreeEigenMatrixMap& jacobianMap) {
        Eigen::Vector3d distance, appliedForce;

        for (size_t i = 0; i < foot.positionPoints.size(); ++i) {
            double deltaXYDerivative = activationXY.evalDerivative(stateVariables(foot.positionPoints[i])(2));

            distance = iDynTree::toEigen(stateVariables(foot.positionPoints[i])) - iDynTree::toEigen(comPosition);
            appliedForce = iDynTree::toEigen(stateVariables(foot.forcePoints[i]));

            jacobianMap.block<3,3>(momentumRange.offset, foot.forcePoints[i].offset).setIdentity();
            jacobianMap.block<3,3>(momentumRange.offset+3, foot.forcePoints[i].offset) = iDynTree::skew(distance);
            jacobianMap.block<3,3>(momentumRange.offset+3, foot.positionPoints[i].offset) = -iDynTree::skew(appliedForce);
            jacobianMap.block<3,3>(momentumRange.offset+3, comPositionRange.offset) += iDynTree::skew(appliedForce);
            jacobianMap.block<2,1>(foot.positionPoints[i].offset, foot.positionPoints[i].offset + 2) = deltaXYDerivative * iDynTree::toEigen(controlVariables(foot.velocityControlPoints[i])).topRows<2>();

//            jacobianMap.block<3,3>(momentumRange.offset+3, basePositionRange.offset) += iDynTree::skew(appliedForce) * iDynTree::toEigen(comjacobianBuffer).leftCols<3>();
//            jacobianMap.block<3,4>(momentumRange.offset+3, baseQuaternionRange.offset) += iDynTree::skew(appliedForce) * iDynTree::toEigen(comjacobianBuffer).block<3,3>(0,3) *
//                    iDynTree::toEigen(iDynTree::Rotation::QuaternionRightTrivializedDerivativeInverse(baseQuaternionNormalized)) *
//                    iDynTree::toEigen(NormalizedQuaternionDerivative(baseQuaternion));
//            jacobianMap.block(momentumRange.offset + 3, jointsPositionRange.offset, 3, jointsPositionRange.size) += iDynTree::skew(appliedForce) * iDynTree::toEigen(comjacobianBuffer).rightCols(jointsPositionRange.size);
        }
    }

    void computeFootRelatedControlJacobian(FootRanges& foot, iDynTree::iDynTreeEigenMatrixMap& jacobianMap) {
        for (size_t i = 0; i < foot.positionPoints.size(); ++i) {
            double deltaXY = activationXY.eval(stateVariables(foot.positionPoints[i])(2));

            jacobianMap.block<3,3>(foot.forcePoints[i].offset, foot.forceControlPoints[i].offset).setIdentity();

            jacobianMap.block<3,3>(foot.positionPoints[i].offset, foot.velocityControlPoints[i].offset).setIdentity();
            jacobianMap.block<2,2>(foot.positionPoints[i].offset, foot.velocityControlPoints[i].offset) *= deltaXY;
        }
    }

    void updateRobotState() {

        robotState = sharedKinDyn->currentState();

        iDynTree::toEigen(basePosition) = iDynTree::toEigen(stateVariables(basePositionRange));
        baseQuaternion = stateVariables(baseQuaternionRange);
        baseQuaternionNormalized = NormalizedQuaternion(baseQuaternion);
        assert(QuaternionBoundsRespected(baseQuaternionNormalized));
        baseRotation.fromQuaternion(baseQuaternionNormalized);
        iDynTree::toEigen(baseQuaternionVelocity) = iDynTree::toEigen(controlVariables(baseQuaternionDerivativeRange));

        robotState.world_T_base.setRotation(baseRotation);
        robotState.world_T_base.setPosition(basePosition);

        robotState.s = stateVariables(jointsPositionRange);

        robotState.s_dot = controlVariables(jointsVelocityRange);

        iDynTree::LinVelocity baseLinVelocity, baseAngVelocity;

        iDynTree::toEigen(baseLinVelocity) = iDynTree::toEigen(controlVariables(baseLinearVelocityRange));
        iDynTree::toEigen(baseAngVelocity) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeInverse(baseQuaternionNormalized)) * iDynTree::toEigen(baseQuaternionVelocity);


        robotState.base_velocity = iDynTree::Twist(baseLinVelocity, baseAngVelocity);
    }

    void setFootRelatedStateSparsity(const FootRanges& foot) {

        for (size_t i = 0; i < foot.positionPoints.size(); ++i) {
            stateSparsity.addIdentityBlock(static_cast<size_t>(momentumRange.offset), static_cast<size_t>(foot.forcePoints[i].offset), 3);
            stateSparsity.addDenseBlock(static_cast<size_t>(momentumRange.offset) + 3, static_cast<size_t>(foot.forcePoints[i].offset), 3, 3);
            stateSparsity.addDenseBlock(static_cast<size_t>(momentumRange.offset) + 3, static_cast<size_t>(foot.positionPoints[i].offset), 3, 3);
            stateSparsity.addDenseBlock(static_cast<size_t>(foot.positionPoints[i].offset), static_cast<size_t>(foot.positionPoints[i].offset) + 2, 2, 1);
        }
    }

    void setFootRelatedControlSparsity(const FootRanges& foot) {
        for (size_t i = 0; i < foot.positionPoints.size(); ++i) {
            controlSparsity.addIdentityBlock(static_cast<size_t>(foot.forcePoints[i].offset), static_cast<size_t>(foot.forceControlPoints[i].offset), 3);
            controlSparsity.addIdentityBlock(static_cast<size_t>(foot.positionPoints[i].offset), static_cast<size_t>(foot.velocityControlPoints[i].offset), 3);
        }
    }

    void setSparsity() {
        stateSparsity.clear();
        controlSparsity.clear();

        setFootRelatedStateSparsity(leftRanges);
        setFootRelatedStateSparsity(rightRanges);
        stateSparsity.addDenseBlock(static_cast<size_t>(momentumRange.offset) + 3, static_cast<size_t>(comPositionRange.offset), 3, 3);
        stateSparsity.addIdentityBlock(static_cast<size_t>(comPositionRange.offset), static_cast<size_t>(momentumRange.offset), 3);
        stateSparsity.addDenseBlock(basePositionRange, baseQuaternionRange);
        stateSparsity.addDenseBlock(baseQuaternionRange, baseQuaternionRange);

        setFootRelatedControlSparsity(leftRanges);
        setFootRelatedControlSparsity(rightRanges);
        controlSparsity.addDenseBlock(basePositionRange, baseLinearVelocityRange);
//        controlSparsity.addDenseBlock(static_cast<size_t>(baseQuaternionRange.offset), static_cast<size_t>(baseVelocityRange.offset) + 3, 4, 3);
        controlSparsity.addIdentityBlock(static_cast<size_t>(baseQuaternionRange.offset), static_cast<size_t>(baseQuaternionDerivativeRange.offset), 4);
        controlSparsity.addIdentityBlock(static_cast<size_t>(jointsPositionRange.offset), static_cast<size_t>(jointsVelocityRange.offset), static_cast<size_t>(jointsPositionRange.size));
    }
};


DynamicalConstraints::DynamicalConstraints(const VariablesLabeller &stateVariables, const VariablesLabeller &controlVariables, std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn, const HyperbolicTangent& planarVelocityActivation)
   : iDynTree::optimalcontrol::DynamicalSystem (stateVariables.size(), controlVariables.size())
   , m_pimpl(std::make_unique<Implementation>())
{
    assert(timelySharedKinDyn);
    assert(timelySharedKinDyn->isValid());
    m_pimpl->timedSharedKinDyn = timelySharedKinDyn;
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;
    m_pimpl->dynamics = m_pimpl->stateVariables;
    m_pimpl->dynamics.zero();

    m_pimpl->totalMass = 0.0;

    const iDynTree::Model & model = timelySharedKinDyn->model();

    for(size_t l=0; l < model.getNrOfLinks(); l++)
    {
        m_pimpl->totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
    }


    m_pimpl->gravityVector.zero();
    m_pimpl->gravityVector(2) = -9.81;

//    m_pimpl->stateJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(stateVariables.size()));
//    m_pimpl->stateJacobianBuffer.zero();

//    m_pimpl->controlJacobianBuffer.resize(static_cast<unsigned int>(stateVariables.size()), static_cast<unsigned int>(controlVariables.size()));
//    m_pimpl->controlJacobianBuffer.zero();

    m_pimpl->activationXY = planarVelocityActivation;

    size_t leftPoints = 0, rightPoints = 0;
    for (auto& label : stateVariables.listOfLabels()) {
        if (label.find("LeftForcePoint") != std::string::npos) {
            leftPoints++;
        }

        if (label.find("RightForcePoint") != std::string::npos) {
             rightPoints++;
        }
    }

    m_pimpl->checkFootVariables("Left",leftPoints, m_pimpl->leftRanges);
    m_pimpl->checkFootVariables("Right", rightPoints, m_pimpl->rightRanges);

    m_pimpl->momentumRange = m_pimpl->stateVariables.getIndexRange("Momentum");
    assert(m_pimpl->momentumRange.isValid());

    m_pimpl->comPositionRange = m_pimpl->stateVariables.getIndexRange("CoMPosition");
    assert(m_pimpl->comPositionRange.isValid());

    m_pimpl->basePositionRange = m_pimpl->stateVariables.getIndexRange("BasePosition");
    assert(m_pimpl->basePositionRange.isValid());

    m_pimpl->baseQuaternionRange = m_pimpl->stateVariables.getIndexRange("BaseQuaternion");
    assert(m_pimpl->baseQuaternionRange.isValid());

    m_pimpl->jointsPositionRange = m_pimpl->stateVariables.getIndexRange("JointsPosition");
    assert(m_pimpl->jointsPositionRange.isValid());

//    m_pimpl->baseVelocityRange = m_pimpl->controlVariables.getIndexRange("BaseVelocity");
//    assert(m_pimpl->baseVelocityRange.isValid());

    m_pimpl->baseLinearVelocityRange = m_pimpl->controlVariables.getIndexRange("BaseLinearVelocity");
    assert(m_pimpl->baseLinearVelocityRange.isValid());

    m_pimpl->baseQuaternionDerivativeRange = m_pimpl->controlVariables.getIndexRange("BaseQuaternionDerivative");
    assert(m_pimpl->baseQuaternionDerivativeRange.isValid());

    m_pimpl->jointsVelocityRange = m_pimpl->controlVariables.getIndexRange("JointsVelocity");
    assert(m_pimpl->jointsVelocityRange.isValid());

    m_pimpl->comjacobianBuffer.resize(3, 6 + static_cast<unsigned int>(m_pimpl->jointsPositionRange.size));
    m_pimpl->comjacobianBuffer.zero();

    m_pimpl->setSparsity();
}

DynamicalConstraints::~DynamicalConstraints()
{

}

bool DynamicalConstraints::dynamics(const iDynTree::VectorDynSize &state, double time, iDynTree::VectorDynSize &stateDynamics)
{
    m_pimpl->stateVariables = state; //this line must remain before those computing the feet related quantities
    m_pimpl->controlVariables = controlInput(); //this line must remain before those computing the feet related quantities

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();

//    m_pimpl->comPosition = m_pimpl->sharedKinDyn->getCenterOfMassPosition(m_pimpl->robotState);

    iDynTree::toEigen(m_pimpl->comPosition) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->comPositionRange));
    iDynTree::toEigen(m_pimpl->dynamics(m_pimpl->momentumRange)) = m_pimpl->totalMass * iDynTree::toEigen(m_pimpl->gravityVector); //this line must remain before those computing the feet related quantities

    m_pimpl->computeFootRelatedDynamics(m_pimpl->leftRanges);

    m_pimpl->computeFootRelatedDynamics(m_pimpl->rightRanges);

    iDynTree::toEigen(m_pimpl->dynamics(m_pimpl->comPositionRange)) = iDynTree::toEigen(m_pimpl->stateVariables(m_pimpl->momentumRange)).topRows<3>()/m_pimpl->totalMass;

    iDynTree::toEigen(m_pimpl->dynamics(m_pimpl->basePositionRange)) = iDynTree::toEigen(m_pimpl->baseRotation) * iDynTree::toEigen(m_pimpl->robotState.base_velocity.getLinearVec3());

//    iDynTree::toEigen(m_pimpl->dynamics(m_pimpl->baseQuaternionRange)) = iDynTree::toEigen(QuaternionLeftTrivializedDerivative(m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(m_pimpl->robotState.base_velocity.getAngularVec3());

    iDynTree::toEigen(m_pimpl->dynamics(m_pimpl->baseQuaternionRange)) = iDynTree::toEigen(m_pimpl->baseQuaternionVelocity);

    iDynTree::toEigen(m_pimpl->dynamics(m_pimpl->jointsPositionRange)) = iDynTree::toEigen(m_pimpl->robotState.s_dot);

    stateDynamics = m_pimpl->dynamics.values();

    return true;
}

bool DynamicalConstraints::dynamicsStateFirstDerivative(const iDynTree::VectorDynSize &state, double time, iDynTree::MatrixDynSize &dynamicsDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = controlInput();
    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();
//    m_pimpl->comPosition = m_pimpl->sharedKinDyn->getCenterOfMassPosition(m_pimpl->robotState);
//    bool ok = m_pimpl->sharedKinDyn->getCenterOfMassJacobian(m_pimpl->robotState, m_pimpl->comjacobianBuffer, iDynTree::FrameVelocityRepresentation::MIXED_REPRESENTATION);

//    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->stateJacobianBuffer);
    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(dynamicsDerivative);

//    jacobianMap.block<3,3>(m_pimpl->momentumRange.offset+3, m_pimpl->basePositionRange.offset).setZero();
//    jacobianMap.block<3,4>(m_pimpl->momentumRange.offset+3, m_pimpl->baseQuaternionRange.offset).setZero();
//    jacobianMap.block(m_pimpl->momentumRange.offset + 3, m_pimpl->jointsPositionRange.offset, 3, m_pimpl->jointsPositionRange.size).setZero();
    jacobianMap.block<3,3>(m_pimpl->momentumRange.offset+3, m_pimpl->comPositionRange.offset).setZero();

    m_pimpl->computeFootRelatedStateJacobian(m_pimpl->leftRanges, jacobianMap);
    m_pimpl->computeFootRelatedStateJacobian(m_pimpl->rightRanges, jacobianMap);

    jacobianMap.block<3,3>(m_pimpl->comPositionRange.offset, m_pimpl->momentumRange.offset).setIdentity();
    jacobianMap.block<3,3>(m_pimpl->comPositionRange.offset, m_pimpl->momentumRange.offset) *= 1.0/m_pimpl->totalMass;


    iDynTree::Matrix4x4 normalizedQuaternionDerivative = NormalizedQuaternionDerivative(m_pimpl->baseQuaternion);
    jacobianMap.block<3,4>(m_pimpl->basePositionRange.offset, m_pimpl->baseQuaternionRange.offset) =
            iDynTree::toEigen(RotatedVectorQuaternionJacobian(m_pimpl->robotState.base_velocity.getLinearVec3(), m_pimpl->baseQuaternionNormalized)) * iDynTree::toEigen(normalizedQuaternionDerivative);

//    jacobianMap.block<4,4>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseQuaternionRange.offset) = iDynTree::toEigen(QuaternionLeftTrivializedDerivativeTimesOmegaJacobian(m_pimpl->robotState.base_velocity.getAngularVec3())) * iDynTree::toEigen(normalizedQuaternionDerivative);

   // dynamicsDerivative = m_pimpl->stateJacobianBuffer;
    return true;
}

bool DynamicalConstraints::dynamicsControlFirstDerivative(const iDynTree::VectorDynSize &state, double time, iDynTree::MatrixDynSize &dynamicsDerivative)
{
    m_pimpl->stateVariables = state;
    m_pimpl->controlVariables = controlInput();

    m_pimpl->sharedKinDyn = m_pimpl->timedSharedKinDyn->get(time);

    m_pimpl->updateRobotState();

//    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(m_pimpl->controlJacobianBuffer);
    iDynTree::iDynTreeEigenMatrixMap jacobianMap = iDynTree::toEigen(dynamicsDerivative);


    m_pimpl->computeFootRelatedControlJacobian(m_pimpl->leftRanges, jacobianMap);
    m_pimpl->computeFootRelatedControlJacobian(m_pimpl->rightRanges, jacobianMap);


    jacobianMap.block<3,3>(m_pimpl->basePositionRange.offset, m_pimpl->baseLinearVelocityRange.offset) = iDynTree::toEigen(m_pimpl->baseRotation);

    jacobianMap.block<4,4>(m_pimpl->baseQuaternionRange.offset, m_pimpl->baseQuaternionDerivativeRange.offset).setIdentity();

    jacobianMap.block(m_pimpl->jointsPositionRange.offset, m_pimpl->jointsVelocityRange.offset, m_pimpl->jointsPositionRange.size, m_pimpl->jointsVelocityRange.size).setIdentity();

//    dynamicsDerivative = m_pimpl->controlJacobianBuffer;
    return true;
}

bool DynamicalConstraints::dynamicsStateFirstDerivativeSparsity(iDynTree::optimalcontrol::SparsityStructure &stateSparsity)
{
    stateSparsity = m_pimpl->stateSparsity;
    return true;
}

bool DynamicalConstraints::dynamicsControlFirstDerivativeSparsity(iDynTree::optimalcontrol::SparsityStructure &controlSparsity)
{
    controlSparsity = m_pimpl->controlSparsity;
    return true;
}
