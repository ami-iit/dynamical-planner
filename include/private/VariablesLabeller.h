/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#ifndef DPLANNER_VARIABLESLABELLER_H
#define DPLANNER_VARIABLESLABELLER_H

#include <iDynTree/Core/Span.h>
#include <iDynTree/Core/Utils.h>
#include <vector>
#include <unordered_map>
#include <string>

class VariablesLabeller {
    std::vector<double> m_fullVector;

    typedef std::unordered_map<std::string, iDynTree::IndexRange> LabelMap;
    LabelMap m_labelMap;


public:
     VariablesLabeller();

     ~VariablesLabeller();

     bool addLabel(const std::string& name, size_t dimension);

     iDynTree::Span<double> values();

     iDynTree::Span<double> values(size_t startIndex, size_t endIndex);

     iDynTree::Span<double> values(const std::string& labelName);

     iDynTree::IndexRange getIndexRange(const std::string& labelName);

};

#endif // DPLANNER_VARIABLESLABELLER_H
