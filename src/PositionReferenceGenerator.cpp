/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/PositionReferenceGenerator.h>

namespace DynamicalPlanner {

class PositionReferenceGeneratorData {

public:

    PositionReferenceGeneratorData(size_t desiredPoints);

    ~PositionReferenceGeneratorData();

    std::vector<PositionWithTimeRange> desiredPositions;

};

class TimeVaryingWeight : public iDynTree::optimalcontrol::TimeVaryingVector {

    std::shared_ptr<PositionReferenceGeneratorData> m_data;
    iDynTree::VectorDynSize m_outputWeight, m_increaseFactors;

public:

    TimeVaryingWeight(std::shared_ptr<PositionReferenceGeneratorData> data, double increaseFactorX,
                      double increaseFactorY, double increaseFactorZ);

    ~TimeVaryingWeight() override;

    const iDynTree::VectorDynSize& get(double time, bool& isValid) override;
};

class PositionReference : public iDynTree::optimalcontrol::TimeVaryingPosition {

    std::shared_ptr<PositionReferenceGeneratorData> m_data;
    iDynTree::Position m_zeroPosition;

public:

    PositionReference(std::shared_ptr<PositionReferenceGeneratorData> data);

    ~PositionReference() override;

    const iDynTree::Position& get(double time, bool& isValid) override;
};

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

class PositionReferenceGenerator::Implementation
{
public:
    std::shared_ptr<TimeVaryingWeight> weightPointer;
    std::shared_ptr<PositionReference> positionPointer;
    std::shared_ptr<PositionReferenceGeneratorData> data;
};

PositionReferenceGenerator::PositionReferenceGenerator()
{
    m_pimpl = std::make_unique<PositionReferenceGenerator::Implementation>();
    m_pimpl->data = std::make_shared<PositionReferenceGeneratorData>(0);
}

PositionReferenceGenerator::PositionReferenceGenerator(unsigned int desiredPoints, double increaseFactorX, double increaseFactorY, double increaseFactorZ) {

    m_pimpl = std::make_unique<PositionReferenceGenerator::Implementation>();
    init(desiredPoints, increaseFactorX, increaseFactorY, increaseFactorZ);
}

PositionReferenceGenerator::~PositionReferenceGenerator() {}

void PositionReferenceGenerator::init(unsigned int desiredPoints, double increaseFactorX, double increaseFactorY, double increaseFactorZ)
{
    m_pimpl->data = std::make_shared<PositionReferenceGeneratorData>(desiredPoints);
    m_pimpl->weightPointer = std::make_shared<TimeVaryingWeight>(m_pimpl->data, increaseFactorX, increaseFactorY, increaseFactorZ);
    m_pimpl->positionPointer = std::make_shared<PositionReference>(m_pimpl->data);
}

std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> PositionReferenceGenerator::timeVaryingWeight() {
    return m_pimpl->weightPointer;
}

std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingPosition> PositionReferenceGenerator::timeVaryingReference() {
    return m_pimpl->positionPointer;
}

PositionWithTimeRange &PositionReferenceGenerator::operator[](size_t index) {
    return at(index);
}

const PositionWithTimeRange &PositionReferenceGenerator::operator[](size_t index) const{
    return at(index);
}

PositionWithTimeRange &PositionReferenceGenerator::at(size_t index)
{
    return m_pimpl->data->desiredPositions[index];

}

const PositionWithTimeRange &PositionReferenceGenerator::at(size_t index) const
{
    return m_pimpl->data->desiredPositions[index];
}

void PositionReferenceGenerator::resize(size_t newSize) {
    PositionWithTimeRange zeroElement;
    zeroElement.desiredPosition.zero();
    m_pimpl->data->desiredPositions.resize(newSize, zeroElement);
}


}
