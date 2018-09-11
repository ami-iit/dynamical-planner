/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_SOLVER_H
#define DPLANNER_SOLVER_H

#include <DynamicalPlanner/Settings.h>
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

};

#endif // DPLANNER_SOLVER_H
