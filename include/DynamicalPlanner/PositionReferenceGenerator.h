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

class PositionReferenceGenerator;

class PositionReferenceGeneratorData {

    friend class PositionReferenceGenerator;

    PositionReferenceGeneratorData(size_t desiredPoints);

public:

    ~PositionReferenceGeneratorData();

    std::vector<PositionWithTimeRange> desiredPositions;

};

class TimeVaryingWeight : public iDynTree::optimalcontrol::TimeVaryingVector {

    friend class PositionReferenceGenerator;

    std::shared_ptr<PositionReferenceGeneratorData> m_data;
    iDynTree::VectorDynSize m_outputWeight, m_increaseFactors;

    TimeVaryingWeight(std::shared_ptr<PositionReferenceGeneratorData> data, double increaseFactorX,
                      double increaseFactorY, double increaseFactorZ);

public:

    ~TimeVaryingWeight() override;

    const iDynTree::VectorDynSize& get(double time, bool& isValid) override;
};

class PositionReference : public iDynTree::optimalcontrol::TimeVaryingPosition {

    friend class PositionReferenceGenerator;
    std::shared_ptr<PositionReferenceGeneratorData> m_data;
    iDynTree::Position m_zeroPosition;

    PositionReference(std::shared_ptr<PositionReferenceGeneratorData> data);

public:

    ~PositionReference() override;

    const iDynTree::Position& get(double time, bool& isValid) override;
};

class PositionReferenceGenerator {

    std::shared_ptr<TimeVaryingWeight> m_weightPointer;
    std::shared_ptr<PositionReference> m_positionPointer;
    std::shared_ptr<PositionReferenceGeneratorData> m_data;

public:

    PositionReferenceGenerator(unsigned int desiredPoints, double increaseFactorX,
                                double increaseFactorY, double increaseFactorZ);

    ~PositionReferenceGenerator();

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> timeVaryingWeight();

    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> timeVaryingReference();

    PositionWithTimeRange& operator[](size_t index);

    const PositionWithTimeRange& operator[](size_t index) const;


    void resize(size_t newSize);
};

}

#endif // DPLANNER_POSITIONREFERENCEGENERATOR_H
