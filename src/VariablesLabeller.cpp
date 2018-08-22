/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <iostream>
#include <private/VariablesLabeller.h>

using namespace DynamicalPlanner::Private;

VariablesLabeller::VariablesLabeller()
{
}

VariablesLabeller::~VariablesLabeller()
{
}

bool VariablesLabeller::addLabel(const std::string &name, size_t dimension)
{
    std::pair<LabelMap::const_iterator, bool> output;
    iDynTree::IndexRange newRange;
    newRange.size = static_cast<long>(dimension);
    newRange.offset = static_cast<long>(m_fullVector.size());
    output = m_labelMap.insert(LabelMap::value_type(name, newRange));

    if (!output.second) {
        std::cerr << "[ERROR][VariablesLabeller::addLabel] The label " << name << " seems to be already existing." << std::endl;
        return false;
    }
    m_fullVector.resize(m_fullVector.size() + dimension, 0.0);
    m_labelsList.push_back(name);
    return true;
}

size_t VariablesLabeller::size() const
{
    return m_fullVector.size();
}

iDynTree::Span<double> VariablesLabeller::values()
{
    return iDynTree::make_span(m_fullVector);
}

iDynTree::Span<double> VariablesLabeller::values(size_t startIndex, size_t endIndex)
{
    return iDynTree::make_span(m_fullVector).subspan(static_cast<long>(startIndex), static_cast<long>(endIndex - startIndex+1));
}

iDynTree::Span<double> VariablesLabeller::operator()(const std::string &labelName)
{
    LabelMap::const_iterator label = m_labelMap.find(labelName);

    long offset = 0;
    long dimension = 0;

    if (label != m_labelMap.cend()) {
        offset = label->second.offset;
        dimension = label->second.size;
    }

    return iDynTree::make_span(m_fullVector).subspan(offset, dimension);
}

//iDynTree::Span<double> VariablesLabeller::values(const std::string &labelName)
//{
//    LabelMap::const_iterator label = m_labelMap.find(labelName);

//    long offset = 0;
//    long dimension = 0;

//    if (label != m_labelMap.cend()) {
//        offset = label->second.offset;
//        dimension = label->second.size;
//    }

//    return iDynTree::make_span(m_fullVector).subspan(offset, dimension);
//}

iDynTree::IndexRange VariablesLabeller::getIndexRange(const std::string &labelName) const
{
    LabelMap::const_iterator label = m_labelMap.find(labelName);

    if (label != m_labelMap.cend()) {
        return label->second;
    }

    return iDynTree::IndexRange::InvalidRange();
}

size_t VariablesLabeller::numberOfLabels() const
{
    return m_labelMap.size();
}

const std::vector<std::string> &VariablesLabeller::listOfLabels() const
{
    return m_labelsList;
}
