/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_COMPILEDEVALUABLE_H
#define LEVI_COMPILEDEVALUABLE_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>

template<typename GenericsMatrix, typename OutputMatrix>
class levi::CompiledEvaluable {
public:

    CompiledEvaluable() {}

    virtual ~CompiledEvaluable() { }

    virtual void evaluate(const std::vector<GenericsMatrix>& generics, OutputMatrix output) = 0;
};

#endif // LEVI_COMPILEDEVALUABLE_H
