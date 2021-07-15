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

}
}

#endif // DPLANNER_UTILITIES_H
