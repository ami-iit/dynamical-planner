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
#include <iDynTree/Core/VectorDynSize.h>
#include <vector>
#include <unordered_map>
#include <string>
namespace DynamicalPlanner {
    namespace Private {
        class VariablesLabeller;
    }
}
class DynamicalPlanner::Private::VariablesLabeller {
    iDynTree::VectorDynSize m_fullVector;
    std::vector<std::string> m_labelsList;

    typedef std::unordered_map<std::string, iDynTree::IndexRange> LabelMap;
    LabelMap m_labelMap;


public:
     VariablesLabeller();

     ~VariablesLabeller();

     bool addLabel(const std::string& name, size_t dimension);

     iDynTree::IndexRange addLabelAndGetIndexRange(const std::string& name, size_t dimension);

     size_t size() const;

     iDynTree::Span<double> values();

     iDynTree::Span<const double> values() const;

     void zero();

     iDynTree::Span<double> operator()(const iDynTree::IndexRange& indexRange);

     iDynTree::Span<const double> operator()(const iDynTree::IndexRange& indexRange) const;

     iDynTree::Span<double> operator()(const std::string& labelName);

     iDynTree::Span<const double> operator()(const std::string& labelName) const;

     double operator()(unsigned int index) const;

     double& operator()(const unsigned int index);

     DynamicalPlanner::Private::VariablesLabeller& operator=(const iDynTree::VectorDynSize& iDynVector);

     iDynTree::IndexRange getIndexRange(const std::string& labelName) const;

     size_t numberOfLabels() const;

     const std::vector<std::string> & listOfLabels() const;

     void clear();

};

#endif // DPLANNER_VARIABLESLABELLER_H
