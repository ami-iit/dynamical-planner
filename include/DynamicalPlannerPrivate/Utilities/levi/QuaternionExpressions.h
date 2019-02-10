/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_QUATERNIONEXPRESSIONS_H
#define DPLANNER_QUATERNIONEXPRESSIONS_H

#include <levi/ForwardDeclarations.h>
#include <string>
#include <memory>

namespace DynamicalPlanner {
    namespace Private {
        class RotationExpression;

        class E_Expression;

        class G_Expression;
    }
}

class DynamicalPlanner::Private::RotationExpression {
    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    RotationExpression(const levi::Variable& quaternionVariable);

    ~RotationExpression();

    levi::Expression& operator->();

    levi::Expression& operator*();
};

class DynamicalPlanner::Private::E_Expression {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    E_Expression(const levi::Variable& quaternionVariable);

    ~E_Expression();

    levi::Expression& operator->();

    levi::Expression& operator*();
};

class DynamicalPlanner::Private::G_Expression {

    class Implementation;
    std::unique_ptr<Implementation> m_pimpl;

public:

    G_Expression(const levi::Variable& quaternionVariable);

    ~G_Expression();

    levi::Expression& operator->();

    levi::Expression& operator*();
};

#endif // DPLANNER_QUATERNIONEXPRESSIONS_H
