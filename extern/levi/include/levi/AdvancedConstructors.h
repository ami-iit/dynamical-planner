/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_ADVANCEDCONSTRUCTORS_H
#define LEVI_ADVANCEDCONSTRUCTORS_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Evaluable.h>
#include <levi/Variable.h>
#include <levi/OperatorsBase.h>
#include <levi/Expression.h>

/**
 * The ConstructorByRows.
 *
 * Constructs an evaluable by stacking the specified rows in the constructor.
 */
template <typename EvaluableT, int rowsNumber>
class levi::ConstructorByRows : public levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, rowsNumber, EvaluableT::cols_at_compile_time>> {

public:

    typedef Eigen::Matrix<typename EvaluableT::value_type, rowsNumber, EvaluableT::cols_at_compile_time> composite_matrix_type;

    typedef levi::Evaluable<composite_matrix_type> composite_evaluable;

private:

    std::vector<levi::ExpressionComponent<EvaluableT>> m_rows;
    std::vector<levi::ExpressionComponent<typename EvaluableT::derivative_evaluable>> m_derivatives;

    template<bool value>
    void eval_row(levi::bool_value<value>, Eigen::Index row);

    void eval_row(levi::bool_value<true>, Eigen::Index row) {
        this->m_evaluationBuffer(row, 0) = m_rows[row].evaluate();
    }

    void eval_row(levi::bool_value<false>, Eigen::Index row) {
        this->m_evaluationBuffer.row(row) = m_rows[row].evaluate();
    }

public:

    ConstructorByRows(const std::vector<levi::ExpressionComponent<EvaluableT>>& rows, std::string name)
        : composite_evaluable(name)
    {
        assert(rows.size() != 0);
        assert((rowsNumber == Eigen::Dynamic) || (rowsNumber == rows.size()));
        Eigen::Index nCols;

        m_rows.push_back(rows.front());
        nCols = m_rows.front().cols();
        this->addDependencies(rows.front());

        for (size_t i = 1; i < rows.size(); ++i) {
            m_rows.push_back(rows[i]);
            assert(m_rows[i].cols() == nCols);
            this->addDependencies(rows[i]);
        }

        this->m_evaluationBuffer.resize(m_rows.size(), nCols);

        m_derivatives.resize(m_rows.size());
    }

    virtual ~ConstructorByRows() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename composite_evaluable::row_type>> row(Eigen::Index row) final {
        return m_rows[row];
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename composite_evaluable::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return m_rows[row](0, col);
    }

    virtual const typename composite_evaluable::matrix_type& evaluate() final {
        for (size_t i = 0; i < m_rows.size(); ++i) {
            eval_row(levi::bool_value<std::is_arithmetic<typename EvaluableT::matrix_type>::value>(), i);
        }
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename composite_evaluable::derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {

        for (size_t i = 0; i < m_rows.size(); ++i) {
            m_derivatives[i] = m_rows[i](0, column).getColumnDerivative(0, variable); //the i-th row of the column derivative corresponds to the (only) column derivative of the element (i, column)
        }

        levi::ExpressionComponent<typename composite_evaluable::derivative_evaluable> derivative;

        derivative = ComposeByRows(m_derivatives, "d(" + this->name() + ")/d" + variable->variableName());

        return derivative;
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        for (auto& expression : m_rows) {
            expression.clearDerivativesCache();
        }
    }

};
template <typename EvaluableT, int rowsNumber>
levi::ConstructorByRows<EvaluableT, rowsNumber>::~ConstructorByRows() { }

/**
 * The ConstructorByCols.
 *
 * Constructs an evaluable by aligning the columns specified in the constructor.
 */
template <typename EvaluableT, int colsNumber>
class levi::ConstructorByCols : public levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, EvaluableT::rows_at_compile_time, colsNumber>> {

public:

    typedef Eigen::Matrix<typename EvaluableT::value_type, EvaluableT::rows_at_compile_time, colsNumber> composite_matrix_type;

    typedef levi::Evaluable<composite_matrix_type> composite_evaluable;

private:

    std::vector<levi::ExpressionComponent<EvaluableT>> m_cols;

    template<bool value>
    void eval_col(levi::bool_value<value>, Eigen::Index col);

    void eval_col(levi::bool_value<true>, Eigen::Index col) {
        this->m_evaluationBuffer(0, col) = m_cols[col].evaluate();
    }

    void eval_col(levi::bool_value<false>, Eigen::Index col) {
        this->m_evaluationBuffer.col(col) = m_cols[col].evaluate();
    }

public:

    ConstructorByCols(const std::vector<levi::ExpressionComponent<EvaluableT>>& cols, std::string name)
        : composite_evaluable(name)
    {
        assert(cols.size() != 0);
        assert((colsNumber == Eigen::Dynamic) || (colsNumber == cols.size()));
        Eigen::Index nRows;

        m_cols.push_back(cols.front());
        nRows = m_cols.front().rows();
        this->addDependencies(cols.front());

        for (size_t i = 1; i < cols.size(); ++i) {
            m_cols.push_back(cols[i]);
            assert(m_cols[i].rows() == nRows);
            this->addDependencies(cols[i]);
        }

        this->m_evaluationBuffer.resize(nRows, m_cols.size());
    }

    virtual ~ConstructorByCols() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename composite_evaluable::col_type>> col(Eigen::Index col) final {
        return m_cols[col];
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename composite_evaluable::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return m_cols[col](row, 0);
    }

    virtual const typename composite_evaluable::matrix_type& evaluate() final {
        for (size_t i = 0; i < m_cols.size(); ++i) {
            eval_col(levi::bool_value<std::is_arithmetic<typename EvaluableT::matrix_type>::value>(), i);
        }
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename composite_evaluable::derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {

        levi::ExpressionComponent<typename composite_evaluable::derivative_evaluable> derivative;

        derivative = m_cols[column].getColumnDerivative(0, variable);

        return derivative;
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        for (auto& expression : m_cols) {
            expression.clearDerivativesCache();
        }
    }

};
template <typename EvaluableT, int colsNumber>
levi::ConstructorByCols<EvaluableT, colsNumber>::~ConstructorByCols(){}

template <typename EvaluableT>
class levi::VariableFromExpressionEvaluable : public levi::EvaluableVariable<typename EvaluableT::col_type> {
    levi::ExpressionComponent<EvaluableT> m_expression;
    bool m_isDependent;

public:

    VariableFromExpressionEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, bool treatAsDependentExpression)
        : levi::EvaluableVariable<typename EvaluableT::col_type>(expression.rows(), expression.name())
          , m_expression(expression)
          , m_isDependent(treatAsDependentExpression)
    {
        static_assert (EvaluableT::cols_at_compile_time == 1 || EvaluableT::cols_at_compile_time == Eigen::Dynamic, "Can't obtain a variable from a matrix evaluable" );
        assert(expression.cols() == 1 && "Can't obtain a variable from a matrix evaluable");
        this->addDependencies(expression);
    }

    virtual ~VariableFromExpressionEvaluable() final;

    virtual const typename EvaluableT::col_type& evaluate() final {
        this->m_evaluationBuffer = m_expression.evaluate();
        return this->m_evaluationBuffer;
    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) final{
        return ((this->variableName() == variable->variableName()) && (this->dimension() == variable->dimension())) || (m_isDependent && m_expression.isDependentFrom(variable));
    }

    virtual std::vector<levi::Registrar*> getDependencies() final {
        std::vector<Registrar*> deps;
        if (m_isDependent) {
            for (auto& dep : this->m_dependencies) {
                deps.push_back(dep.get());
            }
        }
        deps.push_back(this);
        return deps;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::col_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                            std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);
        if ((this->variableName() == variable->variableName()) && (this->dimension() == variable->dimension())) {
            return this->m_identityDerivative;
        } else if (m_isDependent && m_expression.isDependentFrom(variable)) {
            return m_expression.getColumnDerivative(column, variable);
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<typename EvaluableT::col_type>::derivative_evaluable::matrix_type>>(this->dimension(), variable->dimension(),
                                                                                                                                                              "d " + this->variableName() + "/(d " + variable->variableName() + ")");
        }
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_expression.clearDerivativesCache();
    }

};
template <typename EvaluableT>
levi::VariableFromExpressionEvaluable<EvaluableT>::~VariableFromExpressionEvaluable() { }

template <typename CompositeEvaluable, typename LeftEvaluable, typename RightEvaluable>
class levi::HorzcatEvaluable : public levi::BinaryOperator<typename CompositeEvaluable::matrix_type, LeftEvaluable, RightEvaluable>
{
public:

    HorzcatEvaluable(const levi::ExpressionComponent<LeftEvaluable>& lhs,
                     const levi::ExpressionComponent<RightEvaluable>& rhs, const std::string& name)
        : levi::BinaryOperator<typename CompositeEvaluable::matrix_type, LeftEvaluable, RightEvaluable>(lhs, rhs, lhs.rows(),
                                                                                                        lhs.cols() + rhs.cols(), name)
    {
        this->m_info->type = levi::EvaluableType::Horzcat;
        this->m_info->lhs = lhs;
        this->m_info->rhs = rhs;
    }

    virtual ~HorzcatEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::row_type>> row(Eigen::Index row) final {
        return levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::row_type>>::Horzcat(this->m_lhs.row(row),
                                                                                                          this->m_rhs.row(row),
                                                                                                          "[" + this->name() + "](" +
                                                                                                              std::to_string(row) + ",:)");
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::col_type>> col(Eigen::Index col) final {
        if (col < this->m_lhs.cols()) {
            return this->m_lhs.col(col);
        } else {
            return this->m_rhs.col(col - this->m_lhs.cols());
        }
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        if (col < this->m_lhs.cols()) {
            return this->m_lhs(row, col);
        } else {
            return this->m_rhs(row, col - this->m_lhs.cols());
        }
    }

    virtual const typename CompositeEvaluable::matrix_type& evaluate() final {
        this->m_evaluationBuffer.leftCols(this->m_lhs.cols()).lazyAssign(this->m_lhs.evaluate());
        this->m_evaluationBuffer.rightCols(this->m_rhs.cols()).lazyAssign(this->m_rhs.evaluate());
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename CompositeEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
        if (column < this->m_lhs.cols()) {
            return this->m_lhs.getColumnDerivative(column, variable);
        } else {
            return this->m_rhs.getColumnDerivative(column - this->m_lhs.cols(), variable);
        }
    }
};
template <typename CompositeEvaluable, typename LeftEvaluable, typename RightEvaluable>
levi::HorzcatEvaluable<CompositeEvaluable, LeftEvaluable, RightEvaluable>::~HorzcatEvaluable() { }

template <typename CompositeEvaluable, typename TopEvaluable, typename BottomEvaluable>
class levi::VertcatEvaluable : public levi::BinaryOperator<typename CompositeEvaluable::matrix_type, TopEvaluable, BottomEvaluable>
{
public:

    VertcatEvaluable(const levi::ExpressionComponent<TopEvaluable>& top,
                     const levi::ExpressionComponent<BottomEvaluable>& bottom, const std::string& name)
        : levi::BinaryOperator<typename CompositeEvaluable::matrix_type, TopEvaluable, BottomEvaluable>(top, bottom, top.rows() + bottom.rows(),
                                                                                                        top.cols(), name)
    {
        this->m_info->type = levi::EvaluableType::Vertcat;
        this->m_info->lhs = top;
        this->m_info->rhs = bottom;
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::row_type>> row(Eigen::Index row) final {

        if (row < this->m_lhs.rows()) {
            return this->m_lhs.row(row);
        } else {
            return this->m_rhs.row(row - this->m_lhs.rows());
        }
    }

    virtual ~VertcatEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::col_type>> col(Eigen::Index col) final {
        return levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::col_type>>::Vertcat(this->m_lhs.col(col),
                                                                                                          this->m_rhs.col(col),
                                                                                                          "[" + this->name() + "](:," +
                                                                                                              std::to_string(col) + ")");
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename CompositeEvaluable::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        if (row < this->m_lhs.rows()) {
            return this->m_lhs(row, col);
        } else {
            return this->m_rhs(row - this->m_lhs.rows(), col);
        }
    }

    virtual const typename CompositeEvaluable::matrix_type& evaluate() final {
        this->m_evaluationBuffer.topRows(this->m_lhs.rows()) = this->m_lhs.evaluate();
        this->m_evaluationBuffer.bottomRows(this->m_rhs.rows()) = this->m_rhs.evaluate();
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename CompositeEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
        return levi::ExpressionComponent<
            typename CompositeEvaluable::derivative_evaluable>::Vertcat(this->m_lhs.getColumnDerivative(column, variable),
                                                                        this->m_rhs.getColumnDerivative(column, variable),
                                                                        "d[" + this->name() + "(:," + std::to_string(column) + ")]/d"
                                                                            + variable->variableName());
    }
};
template <typename CompositeEvaluable, typename LeftEvaluable, typename RightEvaluable>
levi::VertcatEvaluable<CompositeEvaluable, LeftEvaluable, RightEvaluable>::~VertcatEvaluable() { }

#endif // LEVI_ADVANCEDCONSTRUCTORS_H
