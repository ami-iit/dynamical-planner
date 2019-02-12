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
    return levi::Expression::ComposeByCols({-quaternionImaginaryPart.col(0),
                                            otherColumns.col(0),
                                            otherColumns.col(1),
                                            otherColumns.col(2)}, "E");
}

levi::Expression DynamicalPlanner::Private::G_Expression(const levi::Variable &normalizedQuaternion)
{
    levi::Expression quaternionImaginaryPart, minusQuaternion, skewQuaternion;
    quaternionImaginaryPart = normalizedQuaternion.block(1,0,3,1);
    minusQuaternion = -quaternionImaginaryPart;
    skewQuaternion = minusQuaternion.skew();
    levi::Expression otherColumns = levi::Identity(3,3) * normalizedQuaternion(0,0) + skewQuaternion;
    return  levi::Expression::ComposeByCols({minusQuaternion.col(0),
                                            otherColumns.col(0),
                                            otherColumns.col(1),
                                            otherColumns.col(2)}, "G");
}
