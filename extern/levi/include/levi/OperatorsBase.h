/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_OPERATORSBASE_H
#define LEVI_OPERATORSBASE_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Expression.h>
#include <levi/VariableBase.h>


template <typename MatrixType, typename EvaluableT>
class levi::UnaryOperator : public levi::Evaluable<MatrixType> {
protected:

    levi::ExpressionComponent<EvaluableT> m_expression;

public:

    UnaryOperator(const levi::ExpressionComponent<EvaluableT>& expression, const std::string& name)
        : levi::Evaluable<MatrixType>(name)
        , m_expression(expression)
    {
        this->addDependencies(expression);
    }

    UnaryOperator(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Evaluable<MatrixType>(rows, cols, name)
        , m_expression(expression)
    {
        this->addDependencies(expression);
    }

    UnaryOperator(const levi::ExpressionComponent<EvaluableT>& expression, const MatrixType& initialValue, const std::string& name)
        : levi::Evaluable<MatrixType>(initialValue, name)
        , m_expression(expression)
    {
        this->addDependencies(expression);
    }

    virtual ~UnaryOperator() override;

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_expression.clearDerivativesCache();
    }

};
template <typename MatrixType, typename EvaluableT>
levi::UnaryOperator<MatrixType, EvaluableT>::~UnaryOperator() { }

template <typename MatrixType, class LeftEvaluable, class RightEvaluable>
class levi::BinaryOperator : public levi::Evaluable<MatrixType> {
protected:

    levi::ExpressionComponent<LeftEvaluable> m_lhs;
    levi::ExpressionComponent<RightEvaluable> m_rhs;

public:

    BinaryOperator(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs, const std::string& name)
        : levi::Evaluable<MatrixType>(name)
        , m_lhs(lhs)
        , m_rhs(rhs)
    {
        this->addDependencies(m_lhs, m_rhs);
    }

    BinaryOperator(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs, Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Evaluable<MatrixType>(rows, cols, name)
        , m_lhs(lhs)
        , m_rhs(rhs)
    {
        this->addDependencies(m_lhs, m_rhs);
    }

    BinaryOperator(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs, const MatrixType& initialValue, const std::string& name)
        : levi::Evaluable<MatrixType>(initialValue, name)
        , m_lhs(lhs)
        , m_rhs(rhs)
    {
        this->addDependencies(m_lhs, m_rhs);
    }

    ~BinaryOperator() override;

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_lhs.clearDerivativesCache();
        m_rhs.clearDerivativesCache();
    }
};
template <typename MatrixType, class LeftEvaluable, class RightEvaluable>
levi::BinaryOperator<MatrixType, LeftEvaluable, RightEvaluable>::~BinaryOperator() { }

#endif // LEVI_OPERATORSBASE_H
