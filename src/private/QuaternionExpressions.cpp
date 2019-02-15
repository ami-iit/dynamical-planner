/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>

levi::Expression DynamicalPlanner::Private::RotationExpression(const levi::Variable &normalizedQuaternion)
{
    levi::Expression quaternionImaginaryPart, skewQuaternion;
    quaternionImaginaryPart = normalizedQuaternion.block(1,0,3,1);
    skewQuaternion = quaternionImaginaryPart.skew();
    levi::Expression twoSkewQuaternion = 2.0 * skewQuaternion;
    return levi::Identity(3,3) + normalizedQuaternion(0,0) * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;
}

levi::Expression DynamicalPlanner::Private::E_Expression(const levi::Variable &normalizedQuaternion)
{
    levi::Expression quaternionImaginaryPart, skewQuaternion;
    quaternionImaginaryPart = normalizedQuaternion.block(1,0,3,1);
    skewQuaternion = quaternionImaginaryPart.skew();
    levi::Expression otherColumns = levi::Identity(3,3) * normalizedQuaternion(0,0) + skewQuaternion;
    return levi::Expression::Horzcat(-quaternionImaginaryPart, otherColumns, "E");
}

levi::Expression DynamicalPlanner::Private::G_Expression(const levi::Variable &normalizedQuaternion)
{
    levi::Expression quaternionImaginaryPart, minusQuaternion, skewQuaternion;
    quaternionImaginaryPart = normalizedQuaternion.block(1,0,3,1);
    minusQuaternion = -quaternionImaginaryPart;
    skewQuaternion = minusQuaternion.skew();
    levi::Expression otherColumns = levi::Identity(3,3) * normalizedQuaternion(0,0) + skewQuaternion;
    return  levi::Expression::Horzcat(minusQuaternion, otherColumns, "G");
}

levi::Expression DynamicalPlanner::Private::NormalizedQuaternion(const levi::Variable &notNormalizedQuaternion)
{
    levi::Expression inputQuaternion = notNormalizedQuaternion;
    return inputQuaternion/(inputQuaternion.transpose() * inputQuaternion).pow(0.5);
}

levi::Expression DynamicalPlanner::Private::BodyTwistFromQuaternionVelocity(const levi::Variable &linearVelocity, const levi::Variable &quaternionVelocity,
                                                                            const levi::Variable &normalizedQuaternion, const std::string& name)
{
    return levi::Expression::Vertcat(linearVelocity, 2.0 * G_Expression(normalizedQuaternion) * quaternionVelocity, name);
}
