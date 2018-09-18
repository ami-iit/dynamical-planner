/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SOLVER_H
#define DPLANNER_SOLVER_H


#include <DynamicalPlanner/Settings.h>
#include <DynamicalPlanner/State.h>
#include <DynamicalPlanner/Control.h>

#include <iDynTree/Optimizer.h>
#include <iDynTree/Integrator.h>

#include <memory>

namespace DynamicalPlanner {
    class Solver;
}

class DynamicalPlanner::Solver{

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    Solver();

    ~Solver();

    bool specifySettings(const Settings& settings);

    bool setInitialState(const State& initialState);

    bool setGuesses(std::shared_ptr<TimeVaryingState> stateGuesses, std::shared_ptr<TimeVaryingControl> controlGuesses);

    bool setOptimizer(std::shared_ptr<iDynTree::optimization::Optimizer> optimizer);

    bool setIntegrator(std::shared_ptr<iDynTree::optimalcontrol::integrators::Integrator> integrationMethod);

    bool solve(std::vector<State>& optimalStates, std::vector<Control>& optimalControls);

    const std::vector<State>& optimalStates() const; //Call these methods also to allocate memory at configuration time

    const std::vector<Control>& optimalControls() const;

};

#endif // DPLANNER_SOLVER_H
