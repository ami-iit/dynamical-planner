/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_RECTANGULARFOOT_H
#define DPLANNER_RECTANGULARFOOT_H

#include <iDynTree/Core/Position.h>
#include <iDynTree/Core/Transform.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <iDynTree/Core/Wrench.h>
#include <iDynTree/Core/GeomVector3.h>
#include <vector>
#include <memory>

namespace DynamicalPlanner {
    class RectangularFoot;
}


class DynamicalPlanner::RectangularFoot {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    RectangularFoot();

    ~RectangularFoot();

    bool setFoot(double xLength, double yLength, const iDynTree::Position& topLeftPointPosition); //The z-direction is supposed to be perpendicular to the foot surface. X is considered forward

    bool getPoints(const iDynTree::Transform& frameTransform, std::vector<iDynTree::Position>& pointsPosition) const;

    void getNormalRatiosFromCoP(double copX, double copY, std::vector<double>& pointsRatios); //The CoP is expressed in frame centered in the foot, i.e. not necessarily the foot frame.

    bool getForces(const iDynTree::Wrench& footWrench, std::vector<iDynTree::Force>& pointsForces); //Friction constraints may not satisfied in some points even if they are satisfied in the footWrench. The footWrench is supposed to be in foot coordinates. The output is in foot coordinates.
};


#endif // DPLANNER_RECTANGULARFOOT_H
