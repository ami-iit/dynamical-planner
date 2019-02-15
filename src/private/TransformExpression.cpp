/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/TransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativePositionExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeQuaternionExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeJacobianExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>

using namespace DynamicalPlanner::Private;

class TransformExpression::Implementation {
public:
    levi::Expression position;
    levi::Expression rotation;
};

TransformExpression::TransformExpression()
    : m_pimpl(std::make_unique<Implementation>())
{ }

TransformExpression::TransformExpression(const levi::Expression &position, const levi::Expression &rotation)
    : m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->position = position;
    m_pimpl->rotation = rotation;

}

TransformExpression::TransformExpression(const TransformExpression &other)
    : m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->position = other.position();
    m_pimpl->rotation = other.rotation();
}

TransformExpression::~TransformExpression()
{ }

levi::Expression TransformExpression::position() const
{
    return m_pimpl->position;
}

levi::Expression TransformExpression::rotation() const {
    return m_pimpl->rotation;
}

levi::Expression TransformExpression::operator*(const levi::Expression &position) const
{
    return m_pimpl->position + m_pimpl->rotation * position;
}

void TransformExpression::operator=(const TransformExpression &other)
{
    m_pimpl->position = other.position();
    m_pimpl->rotation = other.rotation();
}

TransformExpression TransformExpression::RelativeTransform(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn, RobotState *robotState, const std::string &baseFrame, const std::string &targetFrame, levi::Variable jointsVariable, levi::ScalarVariable timeVariable)
{

    levi::Expression jacobian = RelativeLeftJacobianExpression(sharedKinDyn, robotState, baseFrame, targetFrame, jointsVariable, timeVariable);
    levi::Expression quaternion = RelativeQuaternionExpression(sharedKinDyn, robotState, baseFrame, targetFrame, jointsVariable, timeVariable, jacobian);
    levi::Expression rotation = RotationExpression(quaternion.asVariable());
    levi::Expression position = RelativePositionExpression(sharedKinDyn, robotState, baseFrame, targetFrame, jointsVariable, timeVariable, jacobian, rotation);

    return TransformExpression(position, rotation);
}

TransformExpression operator*(const TransformExpression &lhs, const TransformExpression &rhs)
{
    return TransformExpression(lhs.position() + lhs.rotation() * rhs.position(), lhs.rotation() * rhs.rotation());
}
