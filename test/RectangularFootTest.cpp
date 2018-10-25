/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/RectangularFoot.h>
#include <iDynTree/Core/TestUtils.h>
#include <cmath>

int main() {

    DynamicalPlanner::RectangularFoot foot;

    double d = 0.08;
    double l = 0.188;

    iDynTree::Position topLeftPosition(0.125,  0.04, 0.0);
    bool ok = foot.setFoot(l, d, topLeftPosition);
    ASSERT_IS_TRUE(ok);

    std::vector<iDynTree::Position> pointsPosition;
    std::vector<iDynTree::Force> pointsForces;
    ok = foot.getPoints(iDynTree::Transform::Identity(), pointsPosition);
    ASSERT_IS_TRUE(ok);

    iDynTree::Wrench appliedWrenchInCenter, appliedWrench;
    iDynTree::Transform frameToCenterTransform;

    iDynTree::Position centerInTopLeftCoordinates(-(l/2.0), -(d/2.0), 0.0);
    iDynTree::Position centerPositionInFrame = topLeftPosition + centerInTopLeftCoordinates;
    frameToCenterTransform.setRotation(iDynTree::Rotation::Identity());
    frameToCenterTransform.setPosition(centerPositionInFrame);

    double t;
    double fx, fy, fz, tx, ty, tz;
    fz = 150;
    double omega = 2;

    for (size_t i = 0 ; i < 1000; ++i) {
        t = 0.01 * i;
        tx = std::min(d * 0.99 * fz/2, std::max(-d * 0.99 * fz/2, t * d * fz/10.0 * std::sin(omega * t)));
        ty = std::min(l * fz/2, std::max(-l * fz/2, t * l * fz/10.0 * std::cos(omega * t)));
        tz = iDynTree::getRandomDouble(-1.0, 1.0);
        fx = iDynTree::getRandomDouble(-10.0, 10.0);
        fy = iDynTree::getRandomDouble(-10.0, 10.0);

        appliedWrenchInCenter(0) = fx;
        appliedWrenchInCenter(1) = fy;
        appliedWrenchInCenter(2) = fz;
        appliedWrenchInCenter(3) = tx;
        appliedWrenchInCenter(4) = ty;
        appliedWrenchInCenter(5) = tz;

        appliedWrench = frameToCenterTransform * appliedWrenchInCenter;

        ok = foot.getForces(appliedWrench, pointsForces);

        iDynTree::Wrench totalWrench, pointWrench;
        totalWrench.zero();
        pointWrench.zero();
        iDynTree::Transform pointTransform;
        pointTransform.setRotation(iDynTree::Rotation::Identity());

        for (size_t f = 0; f < 4; ++f) {
            pointTransform.setPosition(pointsPosition[f]);
            pointWrench.setLinearVec3(pointsForces[f]);
            totalWrench =totalWrench + pointTransform * pointWrench;
            ASSERT_IS_TRUE(pointsForces[f](2) >= 0);
        }

        ASSERT_EQUAL_SPATIAL_FORCE_TOL(appliedWrench, totalWrench, 1e-5);

    }

    return EXIT_SUCCESS;
}
