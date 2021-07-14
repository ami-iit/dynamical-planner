/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/PositionReferenceGenerator.h>

namespace DynamicalPlanner {

PositionReferenceGeneratorData::PositionReferenceGeneratorData(size_t desiredPoints) {

    desiredPositions.resize(desiredPoints);

    for (auto& el : desiredPositions) {
        el.desiredPosition.zero();
    }

}

PositionReferenceGeneratorData::~PositionReferenceGeneratorData() {}

TimeVaryingWeight::TimeVaryingWeight(std::shared_ptr<PositionReferenceGeneratorData> data, double increaseFactorX, double increaseFactorY, double increaseFactorZ) {
    m_data = data;
    m_increaseFactors.resize(3);
    m_increaseFactors(0) = increaseFactorX;
    m_increaseFactors(1) = increaseFactorY;
    m_increaseFactors(2) = increaseFactorZ;

    m_outputWeight.resize(3);
    m_outputWeight.zero();
}

DynamicalPlanner::TimeVaryingWeight::~TimeVaryingWeight(){}

const iDynTree::VectorDynSize &TimeVaryingWeight::get(double time, bool &isValid) {

    isValid = true;
    std::vector<PositionWithTimeRange>::reverse_iterator activeElement;
    activeElement = std::find_if(m_data->desiredPositions.rbegin(),
                                 m_data->desiredPositions.rend(),
                                 [time](const PositionWithTimeRange& a) -> bool { return a.activeRange.isInRange(time); }); //find the last element in the vector with init time lower than the specified time

    if (activeElement == m_data->desiredPositions.rend()) {
        m_outputWeight.zero();
        return m_outputWeight;
    }

    double increaseAmount = (time - activeElement->activeRange.initTime())/(activeElement->activeRange.endTime() - activeElement->activeRange.initTime());

    for (unsigned int i = 0; i < 3; ++i) {
        m_outputWeight(i) = (m_increaseFactors(i) * increaseAmount + 1.0) * (m_increaseFactors(i) * increaseAmount + 1.0);
    }

    return m_outputWeight;
}

PositionReference::PositionReference(std::shared_ptr<PositionReferenceGeneratorData> data) {
    m_data = data;
}

const iDynTree::Position &PositionReference::get(double time, bool &isValid) {

    isValid = true;
    if (!m_data->desiredPositions.size()) {
        return m_zeroPosition;
    }

    std::vector<PositionWithTimeRange>::reverse_iterator activeElement;
    activeElement = std::find_if(m_data->desiredPositions.rbegin(),
                                 m_data->desiredPositions.rend(),
                                 [time](const PositionWithTimeRange& a) -> bool { return a.activeRange.isInRange(time); }); //find the last element in the vector with init time lower than the specified time

    if (activeElement == m_data->desiredPositions.rend()) {
        return m_zeroPosition;
    }

    return activeElement->desiredPosition;
}

DynamicalPlanner::PositionReference::~PositionReference(){}

PositionReferenceGenerator::PositionReferenceGenerator(unsigned int desiredPoints, double increaseFactorX, double increaseFactorY, double increaseFactorZ) {
    m_data.reset(new PositionReferenceGeneratorData(desiredPoints));
    m_weightPointer.reset(new TimeVaryingWeight(m_data, increaseFactorX, increaseFactorY, increaseFactorZ));
    m_positionPointer.reset(new PositionReference(m_data));
}

PositionReferenceGenerator::~PositionReferenceGenerator() {}

std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> PositionReferenceGenerator::timeVaryingWeight() {
    return m_weightPointer;
}

std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> PositionReferenceGenerator::timeVaryingReference() {
    return m_positionPointer;
}

PositionWithTimeRange &PositionReferenceGenerator::operator[](size_t index) {
    return m_data->desiredPositions[index];
}

void PositionReferenceGenerator::resize(size_t newSize) {
    PositionWithTimeRange zeroElement;
    zeroElement.desiredPosition.zero();
    m_data->desiredPositions.resize(newSize, zeroElement);
}

const PositionWithTimeRange &PositionReferenceGenerator::operator[](size_t index) const {
    return m_data->desiredPositions[index];
}


}
