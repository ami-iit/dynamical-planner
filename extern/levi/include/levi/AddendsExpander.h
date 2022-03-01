/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_ADDENDSEXPANDER_H
#define LEVI_ADDENDSEXPANDER_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Evaluable.h>
#include <levi/TypeDetector.h>

template<typename EvaluableT>
class levi::AddendsExpander {

    using Type = levi::EvaluableType;
    using AddendExpression = levi::ExpressionComponent<typename EvaluableT::EvaluableInfo::operands_evaluable>;

    struct Addend {
        int sign = +1;
        AddendExpression lhs, rhs, initial;
    };

    std::vector<Addend> m_addends;
    std::vector<size_t> m_similarAddends;
    std::vector<size_t> m_differentAddends;


    AddendExpression m_newLhs, m_newRhs;

    void expandAddends(const AddendExpression& addend, int sign = +1) {
        assert(addend.isValidExpression());
        Type type;
        type = addend.info().type;

        if (type == Type::Sum) {
            expandAddends(addend.info().lhs, sign);
            expandAddends(addend.info().rhs, sign);
            return;
        }

        if (type == Type::Subtraction) {
            expandAddends(addend.info().lhs, sign);
            expandAddends(addend.info().rhs, -1 * sign);
            return;
        }

        if (type == Type::InvertedSign) {
            expandAddends(addend.info().lhs, -1 * sign);
            return;
        }

        Addend newAddend;


        if (type == Type::Product) {
            newAddend.lhs = addend.info().lhs;
            AddendExpression rhs = addend.info().rhs;
            while (rhs.info().type == Type::Product) {
                newAddend.lhs = newAddend.lhs * rhs.info().lhs;
                rhs = rhs.info().rhs;
            }
            newAddend.rhs = rhs; //extract the rightmost factor
            newAddend.initial = addend;
            newAddend.sign = sign;

            m_addends.push_back(newAddend);
            return;
        }

        newAddend.lhs = levi::ExpressionComponent<levi::IdentityEvaluable<typename EvaluableT::EvaluableInfo::operands_evaluable::matrix_type>>(1,1);
        newAddend.rhs = addend;
        newAddend.initial = addend;
        newAddend.sign = sign;

        m_addends.push_back(newAddend);
    }

    std::pair<Eigen::Index, Eigen::Index> maximumSizeLhs() {
        Eigen::Index maximumRows = m_addends[0].lhs.rows(), maximumCols = m_addends[0].lhs.cols();

        for (size_t i = 0; i < m_similarAddends.size(); ++i) {
            if (m_addends[m_similarAddends[i]].lhs.rows() > maximumRows) {
                maximumRows = m_addends[m_similarAddends[i]].lhs.rows();
            }

            if (m_addends[m_similarAddends[i]].lhs.cols() > maximumCols) {
                maximumCols = m_addends[m_similarAddends[i]].lhs.cols();
            }
        }

        std::pair<Eigen::Index, Eigen::Index> output;
        output.first = maximumRows;
        output.second = maximumCols;

        return output;
    }

    AddendExpression getLhs(size_t index, std::pair<Eigen::Index, Eigen::Index> requiredDimensions) {

        assert(((m_addends[index].lhs.rows() == requiredDimensions.first) && (m_addends[index].lhs.cols() == requiredDimensions.second)) ||
               ((m_addends[index].lhs.rows() == 1) && (m_addends[index].lhs.cols() == 1) && (requiredDimensions.first == requiredDimensions.second)));

        if ((m_addends[index].lhs.rows() == requiredDimensions.first) && (m_addends[index].lhs.cols() == requiredDimensions.second)) {
            return m_addends[index].lhs;
        }

        AddendExpression identity = levi::ExpressionComponent<levi::IdentityEvaluable<typename EvaluableT::EvaluableInfo::operands_evaluable::matrix_type>>(requiredDimensions.first,
                                                                                                                                                            requiredDimensions.second);

        return m_addends[index].lhs * identity;
    }

public:

    AddendsExpander(const AddendExpression& lhs, const AddendExpression& rhs, int rhsSign = +1) {
        expandAddends(lhs);
        expandAddends(rhs, rhsSign);

        for (size_t i = 1; i < m_addends.size(); ++i) {
            if (m_addends[i].rhs == m_addends[0].rhs) { //it will be searched only for similars of the leftmost addend, since other similar addends will be searched recursively when reassembling the newRhs here below
                m_similarAddends.push_back(i);
            } else {
                m_differentAddends.push_back(i);
            }
        }

        if (m_similarAddends.size()) {

            std::pair<Eigen::Index, Eigen::Index> lhsSize = maximumSizeLhs();

            m_newLhs = (m_addends[0].sign > 0) ? getLhs(0, lhsSize) : -getLhs(0, lhsSize);

            for (size_t i = 0; i < m_similarAddends.size(); ++i) {
                m_newLhs = (m_addends[m_similarAddends[i]].sign > 0) ? m_newLhs + getLhs(m_similarAddends[i], lhsSize) : m_newLhs - getLhs(m_similarAddends[i], lhsSize);
            }

            m_newLhs = m_newLhs * m_addends[0].rhs;

            if (m_differentAddends.size()) {
                m_newRhs = (m_addends[m_differentAddends[0]].sign > 0) ? m_addends[m_differentAddends[0]].initial : -m_addends[m_differentAddends[0]].initial;
                for (size_t i = 1; i < m_differentAddends.size(); ++i) {
                    m_newRhs = (m_addends[m_differentAddends[i]].sign > 0) ? m_newRhs + m_addends[m_differentAddends[i]].initial : m_newRhs - m_addends[m_differentAddends[i]].initial;
                }
            }
        }

    }

    const AddendExpression& lhs() const {
        return m_newLhs;
    }

    const AddendExpression& rhs() const {
        return m_newRhs;
    }


};

#endif // LEVI_ADDENDSEXPANDER_H
