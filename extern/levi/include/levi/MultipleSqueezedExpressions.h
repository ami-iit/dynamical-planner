/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_MULTIPLESQUEEZEDEXPRESSIONS_H
#define LEVI_MULTIPLESQUEEZEDEXPRESSIONS_H

#include <levi/ForwardDeclarations.h>
#include <levi/HelpersForwardDeclarations.h>
#include <levi/Expression.h>
#include <levi/TypeDetector.h>
#include <levi/TreeExpander.h>

template <typename EvaluableT>
class levi::MultipleSqueezedExpressions {
public:

    typedef std::unordered_map<std::string, typename EvaluableT::matrix_type> OutputType;
    typedef std::unordered_map<std::string, levi::ExpressionComponent<EvaluableT>> InputType;

private:

    using SqueezedMatrix = typename levi::TreeComponent<EvaluableT>::SqueezedMatrix;

    using Type = levi::EvaluableType;

    std::vector<size_t> m_finalExpressionIndices;
    std::unordered_map<size_t, std::string> m_indicesToNameMap;
    std::vector<levi::TreeComponent<EvaluableT>> m_expandedExpression;
    std::vector<size_t> m_generics;
    OutputType m_results;


public:

    MultipleSqueezedExpressions(const InputType& elements) {

        for (auto& element: elements) {
            m_indicesToNameMap[m_finalExpressionIndices.size()] = element.first;
            m_finalExpressionIndices.emplace_back(levi::expandTree(element.second, true, m_expandedExpression, m_generics));
            m_results[element.first] = SqueezedMatrix::Zero(element.second.rows(), element.second.cols());
        }
    }

    ~MultipleSqueezedExpressions() { }

    const OutputType& evaluate() {

        for (size_t generic : m_generics) {
            m_expandedExpression[generic].buffer = m_expandedExpression[generic].partialExpression.evaluate(); //first evaluate generics
        }

        levi::EvaluableType type;

        for(typename std::vector<levi::TreeComponent<EvaluableT>>::iterator i = m_expandedExpression.begin();
             i != m_expandedExpression.end(); ++i) {
            type = i->type;

            if (type == Type::Sum) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer + m_expandedExpression[i->rhsIndex].buffer); //this is why I need to evaluate the expanded expression in reverse order
            } else if (type == Type::Subtraction) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer - m_expandedExpression[i->rhsIndex].buffer);
            } else if (type == Type::Product) {

                if (m_expandedExpression[i->lhsIndex].cols() != m_expandedExpression[i->rhsIndex].rows()) {
                    if (m_expandedExpression[i->lhsIndex].rows() == 1 && m_expandedExpression[i->lhsIndex].cols() == 1) {
                        i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer(0,0) * m_expandedExpression[i->rhsIndex].buffer);
                    } else {
                        i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer * m_expandedExpression[i->rhsIndex].buffer(0,0));
                    }
                } else {
                    i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer * m_expandedExpression[i->rhsIndex].buffer);
                }

            } else if (type == Type::Division) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer / m_expandedExpression[i->rhsIndex].buffer(0,0));
            } else if (type == Type::InvertedSign) {
                i->buffer.lazyAssign(-m_expandedExpression[i->lhsIndex].buffer);
            } else if (type == Type::Pow) {
                i->buffer(0,0) = std::pow(m_expandedExpression[i->lhsIndex].buffer(0,0), i->exponent);
            } else if (type == Type::Transpose) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer.transpose());
            } else if (type == Type::Row) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer.row(i->block.startRow));
            } else if (type == Type::Column) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer.col(i->block.startCol));
            } else if (type == Type::Element) {
                i->buffer(0,0) = m_expandedExpression[i->lhsIndex].buffer(i->block.startRow, i->block.startCol);
            } else if (type == Type::Block) {
                i->buffer.lazyAssign(m_expandedExpression[i->lhsIndex].buffer.block(i->block.startRow, i->block.startCol, i->block.rows, i->block.cols));
            } else if (type == Type::Horzcat) {
                i->buffer.leftCols(m_expandedExpression[i->lhsIndex].cols()) = m_expandedExpression[i->lhsIndex].buffer;
                i->buffer.rightCols(m_expandedExpression[i->rhsIndex].cols()) = m_expandedExpression[i->rhsIndex].buffer;
            } else if (type == Type::Vertcat) {
                i->buffer.topRows(m_expandedExpression[i->lhsIndex].rows()) = m_expandedExpression[i->lhsIndex].buffer;
                i->buffer.bottomRows(m_expandedExpression[i->rhsIndex].rows()) = m_expandedExpression[i->rhsIndex].buffer;
            }
        }

        for (size_t i = 0; i < m_finalExpressionIndices.size(); ++i) {
            m_results[m_indicesToNameMap[i]] = m_expandedExpression[m_finalExpressionIndices[i]].buffer;
        }

        return m_results;

    }
};

#endif // LEVI_MULTIPLESQUEEZEDEXPRESSIONS_H
