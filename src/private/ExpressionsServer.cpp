/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/CoMInBaseExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeVelocityExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativePositionExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeQuaternionExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeJacobianExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionErrorExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AbsoluteVelocityExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/MomentumInBaseExpression.h>
#include <initializer_list>
#include <iDynTree/Core/EigenHelpers.h>
#include <cassert>
#include <unordered_map>

using namespace DynamicalPlanner::Private;

using ExpressionMap = std::unordered_map<std::string, levi::Expression>;
using TransformsMap = std::unordered_map<std::string, TransformExpression>;

class ExpressionsServer::Implementation {
public:
    std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn;
    SharedKinDynComputationsPointer kinDyn;
    levi::ScalarVariable time = levi::ScalarVariable("t");
    levi::Variable quaternion = levi::Variable(4, "baseQuaternion");
    levi::Expression quaternionNormalized;
    levi::Expression baseRotation;
    levi::Variable basePositionExpr = levi::Variable(3, "aPb");
    levi::Variable baseLinearVelocity = levi::Variable(3, "baseLinVel");
    levi::Variable baseQuaternionVelocity = levi::Variable(4, "baseQuatVel");
    levi::Variable s, s_dot;

    levi::Expression baseTwist;
    TransformExpression worldToBase;
    levi::Expression comInBase;
    levi::Expression crbi;
    ExpressionMap adjointMap, adjointWrenchMap, velocitiesMap, relativePositionsMap,
        relativeQuaternionsMap, relativeRotationsMap, relativeJacobiansMap, quaternionsErrorsMap, motionSubspacesMap, motionSubspacesMatrixMap,
        motionSubspacesWrenchMap, adjointDerivativeMap, adjointWrenchDerivativeMap, absoluteVelocitiesMap,
        absoluteVelocitiesDerivativeMap, linkInertiaMap, linkInertiaInBaseMap, momentumDoubleDerivativeMap, comHessianMap;
    TransformsMap transformsMap;
    RobotState robotState;

    bool first;

    void clearDerivatives(ExpressionMap& map) {
        for (ExpressionMap::iterator it = map.begin(); it != map.end(); ++it) {
            it->second.clearDerivativesCache();
        }
    }
};



ExpressionsServer::ExpressionsServer(std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn)
    : m_pimpl(std::make_unique<Implementation>())
{
    assert(timelySharedKinDyn);
    m_pimpl->timelySharedKinDyn = timelySharedKinDyn;

    m_pimpl->quaternionNormalized = m_pimpl->quaternion/(m_pimpl->quaternion.transpose() * m_pimpl->quaternion).pow(0.5);
    levi::Expression quaternionReal = m_pimpl->quaternionNormalized(0,0);
    levi::Expression quaternionImaginary = m_pimpl->quaternionNormalized.block(1,0,3,1);
    levi::Expression skewQuaternion = quaternionImaginary.skew();
    levi::Expression twoSkewQuaternion = 2.0 * skewQuaternion;
    m_pimpl->baseRotation = levi::Identity(3,3) + quaternionReal * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;

    m_pimpl->time = 0.0;
    m_pimpl->s = levi::Variable(timelySharedKinDyn->model().getNrOfDOFs(), "s");
    m_pimpl->s_dot = levi::Variable(timelySharedKinDyn->model().getNrOfDOFs(), "s_dot");

    levi::Expression minusQuaternion = -quaternionImaginary;
    levi::Expression angularVelocity = 2.0 * (minusQuaternion * m_pimpl->baseQuaternionVelocity(0,0) + (quaternionReal(0,0) * levi::Identity(3,3) - minusQuaternion.skew()) * m_pimpl->baseQuaternionVelocity.block(1,0,3,1));

    m_pimpl->baseTwist = levi::Expression::Vertcat(m_pimpl->baseLinearVelocity, angularVelocity, "baseTwist");

//    m_pimpl->baseTwist = BodyTwistFromQuaternionVelocity(m_pimpl->baseLinearVelocity, m_pimpl->baseQuaternionVelocity,
//                                                         m_pimpl->quaternionNormalized, "baseTwist");

    m_pimpl->worldToBase = TransformExpression(m_pimpl->basePositionExpr, m_pimpl->baseRotation);
    m_pimpl->comInBase = CoMInBaseExpression(this);

    m_pimpl->crbi = levi::Null(6,6);

    for (iDynTree::LinkIndex l = 0; l < static_cast<int>(model().getNrOfLinks()); ++l) {
        m_pimpl->crbi = m_pimpl->crbi + linkInertiaInBase(l) * adjointTransform(model().getLinkName(l), getFloatingBase());
    }

    m_pimpl->first = true;

}

ExpressionsServer::~ExpressionsServer()
{
    //Some Expressions may have cached some derivatives which point to themselves. By clearing the caches we make sure
    //that all the expressions will be deleted
    m_pimpl->clearDerivatives(m_pimpl->adjointMap);
    m_pimpl->clearDerivatives(m_pimpl->adjointWrenchMap);
    m_pimpl->clearDerivatives(m_pimpl->adjointDerivativeMap);
    m_pimpl->clearDerivatives(m_pimpl->adjointWrenchDerivativeMap);
    m_pimpl->clearDerivatives(m_pimpl->velocitiesMap);
    m_pimpl->clearDerivatives(m_pimpl->relativePositionsMap);
    m_pimpl->clearDerivatives(m_pimpl->relativeQuaternionsMap);
    m_pimpl->clearDerivatives(m_pimpl->relativeRotationsMap);
    m_pimpl->clearDerivatives(m_pimpl->relativeJacobiansMap);
    m_pimpl->clearDerivatives(m_pimpl->quaternionsErrorsMap);
    m_pimpl->clearDerivatives(m_pimpl->absoluteVelocitiesMap);
    m_pimpl->clearDerivatives(m_pimpl->absoluteVelocitiesDerivativeMap);
    m_pimpl->clearDerivatives(m_pimpl->linkInertiaInBaseMap);
    m_pimpl->clearDerivatives(m_pimpl->momentumDoubleDerivativeMap);
    m_pimpl->clearDerivatives(m_pimpl->comHessianMap);
    m_pimpl->comInBase.clearDerivativesCache();
    m_pimpl->crbi.clearDerivativesCache();
}

bool ExpressionsServer::updateRobotState(double time, const RobotState &currentState)
{
    m_pimpl->time = time;
    m_pimpl->kinDyn = m_pimpl->timelySharedKinDyn->get(time);

    if (m_pimpl->first || !(m_pimpl->kinDyn->sameState(currentState))) {
        if (!(m_pimpl->kinDyn->updateRobotState(currentState))){
            return false;
        }

        m_pimpl->quaternion = iDynTree::toEigen(currentState.base_quaternion);
        m_pimpl->basePositionExpr = iDynTree::toEigen(currentState.base_position);
        m_pimpl->s = iDynTree::toEigen(currentState.s);
        m_pimpl->s_dot = iDynTree::toEigen(currentState.s_dot);
        m_pimpl->baseLinearVelocity = iDynTree::toEigen(currentState.base_linearVelocity);
        m_pimpl->baseQuaternionVelocity = iDynTree::toEigen(currentState.base_quaternionVelocity);

        m_pimpl->robotState = m_pimpl->kinDyn->currentState();

        m_pimpl->first = false;
    }

    return true;
}

bool ExpressionsServer::updateRobotState(double time)
{
    m_pimpl->time = time;
    m_pimpl->kinDyn = m_pimpl->timelySharedKinDyn->get(time);

    if (m_pimpl->first || !(m_pimpl->kinDyn->sameState(m_pimpl->robotState))) {

        m_pimpl->quaternion = iDynTree::toEigen(m_pimpl->kinDyn->currentState().base_quaternion);
        m_pimpl->basePositionExpr = iDynTree::toEigen(m_pimpl->kinDyn->currentState().base_position);
        m_pimpl->s = iDynTree::toEigen(m_pimpl->kinDyn->currentState().s);
        m_pimpl->s_dot = iDynTree::toEigen(m_pimpl->kinDyn->currentState().s_dot);
        m_pimpl->baseLinearVelocity = iDynTree::toEigen(m_pimpl->kinDyn->currentState().base_linearVelocity);
        m_pimpl->baseQuaternionVelocity = iDynTree::toEigen(m_pimpl->kinDyn->currentState().base_quaternionVelocity);

        m_pimpl->robotState = m_pimpl->kinDyn->currentState();

        m_pimpl->first = false;
    }

    return true;
}

const RobotState &ExpressionsServer::currentState() const
{
    return m_pimpl->robotState;
}

SharedKinDynComputationsPointer ExpressionsServer::currentKinDyn()
{
    return m_pimpl->timelySharedKinDyn->get(m_pimpl->time.evaluate());
}

const iDynTree::Model &ExpressionsServer::model() const
{
    return m_pimpl->timelySharedKinDyn->model();
}

std::string ExpressionsServer::getFloatingBase() const
{
    return m_pimpl->timelySharedKinDyn->getFloatingBase();
}

levi::Expression ExpressionsServer::baseRotation()
{
    return (m_pimpl->baseRotation);
}

levi::Expression ExpressionsServer::normalizedBaseQuaternion()
{
    return (m_pimpl->quaternionNormalized);
}

levi::Variable ExpressionsServer::baseQuaternion()
{
    return (m_pimpl->quaternion);
}

levi::Variable ExpressionsServer::basePosition()
{
    return (m_pimpl->basePositionExpr);
}

levi::Variable ExpressionsServer::baseLinearVelocity()
{
    return (m_pimpl->baseLinearVelocity);
}

levi::Variable ExpressionsServer::baseQuaternionVelocity()
{
    return (m_pimpl->baseQuaternionVelocity);
}

levi::Expression ExpressionsServer::baseTwist()
{
    return (m_pimpl->baseTwist);
}

levi::Variable ExpressionsServer::jointsPosition()
{
    return (m_pimpl->s);
}

levi::Variable ExpressionsServer::jointsVelocity()
{
    return (m_pimpl->s_dot);
}

TransformExpression ExpressionsServer::worldToBase()
{
    return (m_pimpl->worldToBase);
}

levi::Expression ExpressionsServer::comInBase()
{
    return (m_pimpl->comInBase);
}

levi::Expression ExpressionsServer::comInBaseHessian(long column)
{
    std::string label = std::to_string(column);

    ExpressionMap::iterator element = m_pimpl->comHessianMap.find(label);

    if (element != m_pimpl->comHessianMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;
        newElement.second = CoMInBaseJointsDoubleDerivative(this, column);
        auto result = m_pimpl->comHessianMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::adjointTransform(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->adjointMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->adjointMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Identity(6,6);
        } else {
            newElement.second = AdjointTransformExpression(this, baseFrame, targetFrame);
        }
        auto result = m_pimpl->adjointMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::adjointTransformJointsDerivative(const std::string &baseFrame, const std::string &targetFrame, long column)
{
    std::string label = baseFrame + targetFrame + "__" + std::to_string(column);

    ExpressionMap::iterator element = m_pimpl->adjointDerivativeMap.find(label);

    if (element != m_pimpl->adjointDerivativeMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Null(6, m_pimpl->s.rows());
        } else {
            newElement.second = AdjointTransformExpressionJointsDerivative(this, baseFrame, targetFrame, column);
        }
        auto result = m_pimpl->adjointDerivativeMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::adjointTransformWrench(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->adjointWrenchMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->adjointWrenchMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Identity(6, 6);
        } else {
            newElement.second = AdjointTransformWrenchExpression(this, baseFrame, targetFrame);
        }
        auto result = m_pimpl->adjointWrenchMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::adjointTransformWrenchJointsDerivative(const std::string &baseFrame, const std::string &targetFrame, long column)
{
    std::string label = baseFrame + targetFrame + "__" + std::to_string(column);

    ExpressionMap::iterator element = m_pimpl->adjointWrenchDerivativeMap.find(label);

    if (element != m_pimpl->adjointWrenchDerivativeMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Null(6, m_pimpl->s.rows());
        } else {
            newElement.second = AdjointTransformWrenchExpressionJointsDerivative(this, baseFrame, targetFrame, column);
        }
        auto result = m_pimpl->adjointWrenchDerivativeMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::relativePosition(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->relativePositionsMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->relativePositionsMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Null(3, 1);
        } else {
            newElement.second = RelativePositionExpression(this, baseFrame, targetFrame);
        }
        auto result = m_pimpl->relativePositionsMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::relativeQuaternion(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->relativeQuaternionsMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->relativeQuaternionsMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Identity(4,4).col(0);
        } else {
            newElement.second = RelativeQuaternionExpression(this, baseFrame, targetFrame);
        }
        auto result = m_pimpl->relativeQuaternionsMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::relativeRotation(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->relativeRotationsMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->relativeRotationsMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Identity(3, 3);
        } else {
            newElement.second = adjointTransform(baseFrame, targetFrame).block(0, 0, 3 ,3);
        }
        auto result = m_pimpl->relativeRotationsMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

TransformExpression ExpressionsServer::relativeTransform(const std::string &baseFrame, const std::string &targetFrame)
{
    TransformsMap::iterator element = m_pimpl->transformsMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->transformsMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, TransformExpression> newElement;
        newElement.first = baseFrame + targetFrame;
        newElement.second = TransformExpression(relativePosition(baseFrame, targetFrame), relativeRotation(baseFrame, targetFrame));
        auto result = m_pimpl->transformsMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::relativeLeftJacobian(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->relativeJacobiansMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->relativeJacobiansMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Null(6, m_pimpl->s.rows());
        } else {
            newElement.second = RelativeLeftJacobianExpression(this, baseFrame, targetFrame);
        }
        auto result = m_pimpl->relativeJacobiansMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::relativeVelocity(const std::string &baseFrame, const std::string &targetFrame)
{
    ExpressionMap::iterator element = m_pimpl->velocitiesMap.find(baseFrame+targetFrame);

    if (element != m_pimpl->velocitiesMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = baseFrame + targetFrame;
        if (baseFrame == targetFrame) {
            newElement.second = levi::Null(6, 1);
        } else {
            newElement.second = RelativeLeftVelocityExpression(this, baseFrame, targetFrame);
        }
        auto result = m_pimpl->velocitiesMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::absoluteVelocity(const std::string &targetFrame, const levi::Variable &baseTwist)
{
    std::string key = targetFrame + baseTwist.name();
    ExpressionMap::iterator element = m_pimpl->absoluteVelocitiesMap.find(key);

    if (element != m_pimpl->absoluteVelocitiesMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = key;
        if (targetFrame == getFloatingBase()) {
            newElement.second = baseTwist;
        } else {
            newElement.second = AbsoluteLeftVelocityExpression(this, baseTwist, targetFrame);
        }
        auto result = m_pimpl->absoluteVelocitiesMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::absoluteVelocity(const std::string &targetFrame)
{
    return absoluteVelocity(targetFrame, baseTwist().asVariable());
}

levi::Expression ExpressionsServer::absoluteVelocityJointsDerivative(const std::string &targetFrame, const levi::Variable &baseTwist)
{
    std::string label = targetFrame + baseTwist.name();

    ExpressionMap::iterator element = m_pimpl->absoluteVelocitiesDerivativeMap.find(label);

    if (element != m_pimpl->absoluteVelocitiesDerivativeMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;
        if (targetFrame == getFloatingBase()) {
            newElement.second = levi::Null(6, m_pimpl->s.rows());
        } else {
            newElement.second = AbsoluteLeftVelocityJointsDerivativeExpression(this, baseTwist, targetFrame);
        }
        auto result = m_pimpl->absoluteVelocitiesDerivativeMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::quaternionError(const std::string &desiredFrame, const levi::Variable &desiredQuaternion)
{
    ExpressionMap::iterator element = m_pimpl->quaternionsErrorsMap.find(desiredFrame);

    if (element != m_pimpl->quaternionsErrorsMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = desiredFrame;
        newElement.second = QuaternionError(desiredFrame, this, desiredQuaternion);
        auto result = m_pimpl->quaternionsErrorsMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::motionSubSpaceVector(iDynTree::JointIndex joint, iDynTree::LinkIndex parentLink, iDynTree::LinkIndex childLink)
{
    std::string label = model().getJointName(joint) + std::to_string(parentLink) + std::to_string(childLink);
    ExpressionMap::iterator element = m_pimpl->motionSubspacesMap.find(label);

    if (element != m_pimpl->motionSubspacesMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;

        Eigen::Matrix<double, 6,1> motionSubSpace = iDynTree::toEigen(model().getJoint(joint)->getMotionSubspaceVector(0,
                                                                                                                        childLink,
                                                                                                                        parentLink));

        newElement.second = levi::Constant(motionSubSpace, "s_" + std::to_string(joint));
        auto result = m_pimpl->motionSubspacesMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::motionSubSpaceAsCrossProduct(iDynTree::JointIndex joint, iDynTree::LinkIndex parentLink, iDynTree::LinkIndex childLink)
{
    std::string label = model().getJointName(joint) + std::to_string(parentLink) + std::to_string(childLink);
    ExpressionMap::iterator element = m_pimpl->motionSubspacesMatrixMap.find(label);

    if (element != m_pimpl->motionSubspacesMatrixMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;

        Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProduct =
            iDynTree::toEigen(model().getJoint(joint)->getMotionSubspaceVector(0,
                                                                               childLink,
                                                                               parentLink).asCrossProductMatrix());

        newElement.second = levi::Constant(motionSubSpaceAsCrossProduct, "s_" + std::to_string(joint) + "x");
        auto result = m_pimpl->motionSubspacesMatrixMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::motionSubSpaceAsCrossProductWrench(iDynTree::JointIndex joint, iDynTree::LinkIndex parentLink, iDynTree::LinkIndex childLink)
{
    std::string label = model().getJointName(joint) + std::to_string(parentLink) + std::to_string(childLink);
    ExpressionMap::iterator element = m_pimpl->motionSubspacesWrenchMap.find(label);

    if (element != m_pimpl->motionSubspacesWrenchMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;

        Eigen::Matrix<double, 6,6> motionSubSpaceAsCrossProductWrench =
            iDynTree::toEigen(model().getJoint(joint)->getMotionSubspaceVector(0,
                                                                               childLink,
                                                                               parentLink).asCrossProductMatrixWrench());

        newElement.second = levi::Constant(motionSubSpaceAsCrossProductWrench, "s_" + std::to_string(joint) + "x*");
        auto result = m_pimpl->motionSubspacesWrenchMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::linkInertia(iDynTree::LinkIndex link)
{
    std::string label = model().getLinkName(link);
    ExpressionMap::iterator element = m_pimpl->linkInertiaMap.find(label);

    if (element != m_pimpl->linkInertiaMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;

        newElement.second = levi::Constant(iDynTree::toEigen(model().getLink(link)->getInertia().asMatrix()),"I_" + model().getLinkName(link));
        auto result = m_pimpl->linkInertiaMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::linkInertiaInBase(iDynTree::LinkIndex link)
{
    std::string label = model().getLinkName(link);
    ExpressionMap::iterator element = m_pimpl->linkInertiaInBaseMap.find(label);

    if (element != m_pimpl->linkInertiaInBaseMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;

        newElement.second = adjointTransformWrench(getFloatingBase(), model().getLinkName(link)) * linkInertia(link);
        auto result = m_pimpl->linkInertiaInBaseMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}

levi::Expression ExpressionsServer::compositeRigidBodyInertia()
{
    return m_pimpl->crbi;
}

levi::Expression ExpressionsServer::momentumInBaseJointsDoubleDerivative(const levi::Variable &baseTwist, long column)
{
    std::string label = std::to_string(column);

    ExpressionMap::iterator element = m_pimpl->momentumDoubleDerivativeMap.find(label);

    if (element != m_pimpl->momentumDoubleDerivativeMap.end()) {
        return (element->second);
    } else {
        std::pair<std::string, levi::Expression> newElement;
        newElement.first = label;
        newElement.second = MomentumInBaseExpressionJointsDoubleDerivativeExpression(this, baseTwist, column);
        auto result = m_pimpl->momentumDoubleDerivativeMap.insert(newElement);
        assert(result.second);
        return (result.first->second);
    }
}
