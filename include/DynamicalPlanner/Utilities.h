/*
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_UTILITIES_H
#define DPLANNER_UTILITIES_H

#include <DynamicalPlanner/Settings.h>
#include <DynamicalPlanner/RectangularFoot.h>
#include <DynamicalPlanner/State.h>
#include <DynamicalPlanner/PositionReferenceGenerator.h>
#include <iDynTree/TimeVaryingObject.h>
#include <iDynTree/Core/VectorFixSize.h>
#include <memory>

namespace DynamicalPlanner {
namespace Utilities {

bool FillDefaultInitialState(const DynamicalPlanner::Settings& inputSettings, const iDynTree::VectorDynSize& desiredJoints,
                             DynamicalPlanner::RectangularFoot &leftFoot, DynamicalPlanner::RectangularFoot &rightFoot,
                             DynamicalPlanner::State& initialState);

bool SetMinContactPointToZero(const DynamicalPlanner::Settings& inputSettings, DynamicalPlanner::State &initialState);

class TranslatingCoMStateGuess : public DynamicalPlanner::TimeVaryingState {
    DynamicalPlanner::State m_state, m_initialState;
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> m_comReference;
public:

    TranslatingCoMStateGuess(std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> comReference, const DynamicalPlanner::State &initialState);

    ~TranslatingCoMStateGuess() override;

    DynamicalPlanner::State &get(double time, bool &isValid) override;
};

class SimpleWalkingStateMachine
{
    std::shared_ptr<PositionReferenceGenerator> m_references;

    iDynTree::Position m_stepIncrement;
    double m_stepDuration;
    double m_horizon;
    double m_forceThreshold;
    double m_positionThreshold;
    double m_minimumDt;
    bool m_verbose{false};

public:

    bool initialize(const iDynTree::Vector3 startingReference, const iDynTree::Vector3& stepIncrement,
                    double stepDuration, double horizon, double minimumDt,
                    double weightIncreaseX, double weightIncreaseY, double weightIncreaseZ,
                    double forceThreshold = 40, double positionErrorThreshold = 5e-3);

    void setVerbose(bool verbose = true);

    bool advance(const DynamicalPlanner::State& currentState, const DynamicalPlanner::State& futureState);

    std::shared_ptr<PositionReferenceGenerator> references();

};

}
}

#endif // DPLANNER_UTILITIES_H
