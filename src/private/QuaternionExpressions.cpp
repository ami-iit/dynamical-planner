/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>

using namespace DynamicalPlanner::Private;

class RotationExpression::Implementation {
public:
    levi::Variable quaternion;

    levi::Expression baseRotation;
};

RotationExpression::RotationExpression(const levi::Variable &quaternionVariable)
    : m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->quaternion = quaternionVariable;
    levi::Expression quaternionNormalized, quaternionImaginaryPart, skewQuaternion;
    quaternionNormalized = m_pimpl->quaternion/(m_pimpl->quaternion.transpose() * m_pimpl->quaternion).pow(0.5);
    quaternionImaginaryPart = quaternionNormalized.block(1,0,3,1);
    skewQuaternion = quaternionImaginaryPart.skew();
    levi::Expression twoSkewQuaternion = 2.0 * skewQuaternion;
    m_pimpl->baseRotation = levi::Identity(3,3) + quaternionNormalized(0,0) * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;
}

RotationExpression::~RotationExpression()
{ }

levi::Expression &RotationExpression::operator->()
{
    return m_pimpl->baseRotation;
}

levi::Expression &RotationExpression::operator*()
{
    return m_pimpl->baseRotation;
}

class E_Expression::Implementation {
public:
    levi::Variable quaternion;

    levi::Expression E;
};

E_Expression::E_Expression(const levi::Variable &quaternionVariable)
    : m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->quaternion = quaternionVariable;
    levi::Expression quaternionNormalized, quaternionImaginaryPart, skewQuaternion;
    quaternionNormalized = m_pimpl->quaternion/(m_pimpl->quaternion.transpose() * m_pimpl->quaternion).pow(0.5);
    quaternionImaginaryPart = quaternionNormalized.block(1,0,3,1);
    skewQuaternion = quaternionImaginaryPart.skew();
    levi::Expression otherColumns = levi::Identity(3,3) * quaternionNormalized(0,0) + skewQuaternion;
    m_pimpl->E = levi::Expression::ComposeByCols({-quaternionImaginaryPart.col(0),
                                                  otherColumns.col(0),
                                                  otherColumns.col(1),
                                                  otherColumns.col(2)}, "E");
}

E_Expression::~E_Expression()
{ }

levi::Expression &E_Expression::operator->()
{
    return m_pimpl->E;
}

levi::Expression &E_Expression::operator*()
{
    return m_pimpl->E;
}


class G_Expression::Implementation {
public:
    levi::Variable quaternion;

    levi::Expression G;
};

G_Expression::G_Expression(const levi::Variable &quaternionVariable)
    : m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->quaternion = quaternionVariable;
    levi::Expression quaternionNormalized, quaternionImaginaryPart, minusQuaternion, skewQuaternion;
    quaternionNormalized = m_pimpl->quaternion/(m_pimpl->quaternion.transpose() * m_pimpl->quaternion).pow(0.5);
    quaternionImaginaryPart = quaternionNormalized.block(1,0,3,1);
    minusQuaternion = -quaternionImaginaryPart;
    skewQuaternion = minusQuaternion.skew();
    levi::Expression otherColumns = levi::Identity(3,3) * quaternionNormalized(0,0) + skewQuaternion;
    m_pimpl->G = levi::Expression::ComposeByCols({minusQuaternion.col(0),
                                                  otherColumns.col(0),
                                                  otherColumns.col(1),
                                                  otherColumns.col(2)}, "G");
}

G_Expression::~G_Expression()
{ }

levi::Expression &G_Expression::operator->()
{
    return m_pimpl->G;
}

levi::Expression &G_Expression::operator*()
{
    return m_pimpl->G;
}
