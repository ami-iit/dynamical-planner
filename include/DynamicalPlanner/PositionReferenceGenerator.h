/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#ifndef DPLANNER_POSITIONREFERENCEGENERATOR_H
#define DPLANNER_POSITIONREFERENCEGENERATOR_H

#include <iDynTree/Core/Position.h>
#include <iDynTree/TimeRange.h>
#include <iDynTree/TimeVaryingObject.h>
#include <memory>
#include <vector>

namespace DynamicalPlanner {

struct PositionWithTimeRange{
    iDynTree::Position desiredPosition;
    iDynTree::optimalcontrol::TimeRange activeRange;
};

class PositionReferenceGenerator {

    class Implementation;

    std::unique_ptr<Implementation> m_pimpl;

public:

    PositionReferenceGenerator();

    PositionReferenceGenerator(unsigned int desiredPoints, double increaseFactorX,
                               double increaseFactorY, double increaseFactorZ);

    ~PositionReferenceGenerator();

    void init(unsigned int desiredPoints, double increaseFactorX,
              double increaseFactorY, double increaseFactorZ);

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> timeVaryingWeight();

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> timeVaryingReference();

    PositionWithTimeRange& operator[](size_t index);

    const PositionWithTimeRange& operator[](size_t index) const;

    PositionWithTimeRange& at(size_t index);

    const PositionWithTimeRange& at(size_t index) const;

    void resize(size_t newSize);
};

}

#endif // DPLANNER_POSITIONREFERENCEGENERATOR_H
