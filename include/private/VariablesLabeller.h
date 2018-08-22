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
namespace DynamicalPlanner {
    namespace Private {
        class VariablesLabeller;
    }
}
class DynamicalPlanner::Private::VariablesLabeller {
    std::vector<double> m_fullVector;
    std::vector<std::string> m_labelsList;

    typedef std::unordered_map<std::string, iDynTree::IndexRange> LabelMap;
    LabelMap m_labelMap;


public:
     VariablesLabeller();

     ~VariablesLabeller();

     bool addLabel(const std::string& name, size_t dimension);

     size_t size() const;

     iDynTree::Span<double> values();

     iDynTree::Span<double> operator()(const iDynTree::IndexRange& indexRange);

     iDynTree::Span<double> operator()(const std::string& labelName);

     iDynTree::IndexRange getIndexRange(const std::string& labelName) const;

     size_t numberOfLabels() const;

     const std::vector<std::string> & listOfLabels() const;

};

#endif // DPLANNER_VARIABLESLABELLER_H
