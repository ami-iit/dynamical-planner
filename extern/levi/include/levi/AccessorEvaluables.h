/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_ACCESSOREVALUABLES_H
#define LEVI_ACCESSOREVALUABLES_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/BasicEvaluables.h>
#include <levi/Expression.h>
#include <levi/OperatorsBase.h>


/**
 * Helper struct for determining the type resulting from a block extraction. Specialization for a matrix.
 */
template<typename Matrix>
struct levi::dynamic_block_return<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> {
    typedef Eigen::Matrix<typename Matrix::value_type, Eigen::Dynamic, Eigen::Dynamic> type;
};

/**
 * Helper struct for determining the type resulting from a block extraction. Specialization for a scalar.
 */
template<typename Scalar>
struct levi::dynamic_block_return<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> {
    typedef Scalar type;
};

/**
 * Helper struct for determining the type resulting from a block extraction. Specialization for a matrix.
 */
template<typename Matrix, int rows, int cols>
struct levi::fixed_block_return<Matrix, rows, cols, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> {
    typedef Eigen::Matrix<typename Matrix::value_type, rows, cols> type;
};

/**
 * Helper struct for determining the type resulting from a block extraction. Specialization for a scalar.
 */
template<typename Scalar, int rows, int cols>
struct levi::fixed_block_return<Scalar, rows, cols, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> {
    typedef Scalar type;
};

/**
 * Helper struct for determining the type resulting from a transposition. Specialization for a matrix.
 */
template<typename EvaluableT>
struct levi::transpose_type<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type> {
    typedef Eigen::Matrix<typename EvaluableT::value_type, EvaluableT::cols_at_compile_time, EvaluableT::rows_at_compile_time> type;
};

/**
 * Helper struct for determining the type resulting from a transposition. Specialization for a scalar.
 */
template<typename EvaluableT>
struct levi::transpose_type<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type> {
    typedef typename EvaluableT::matrix_type type;
};



/**
 * @brief The RowEvaluable. To retrieve a specified row from an evaluable. Specialization for matrix valued evaluables.
 */
template <typename EvaluableT>
class levi::RowEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename EvaluableT::row_type, EvaluableT>
{
    Eigen::Index m_row;

public:

    RowEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index row)
        : levi::UnaryOperator<typename EvaluableT::row_type, EvaluableT>(expression, 1, expression.cols(), "(" + expression.name() + ")(" + std::to_string(row) + ",:)")
        , m_row(row)
    {
        this->m_info->type = levi::EvaluableType::Row;
        this->m_info->lhs = expression;
        this->m_info->block.startRow = row;
    }

    virtual ~RowEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<typename EvaluableT::row_type>::col_type>> col(Eigen::Index col) final {
        return this->m_expression(m_row, col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<typename EvaluableT::row_type>::value_type>> element(Eigen::Index row, Eigen::Index col) {
        levi::unused(row);
        return this->m_expression(m_row, col);
    }

    virtual const typename EvaluableT::row_type& evaluate() final {
        this->m_evaluationBuffer.lazyAssign(this->m_expression.evaluate().row(m_row));

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::row_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                            std::shared_ptr<levi::VariableBase> variable) final {
        levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::row_type>::derivative_evaluable> newDerivative;

        newDerivative = levi::ExpressionComponent<levi::RowEvaluable<typename EvaluableT::derivative_evaluable>>(this->m_expression.getColumnDerivative(column, variable), m_row);

        return newDerivative;
    }

};
template <typename EvaluableT>
levi::RowEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>::~RowEvaluable(){ }


/**
 * @brief The RowEvaluable. To retrieve a specified row from an evaluable. Specialization for scalar valued evaluables.
 */
template <typename EvaluableT>
class levi::RowEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename EvaluableT::row_type, EvaluableT>
{

public:

    RowEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index row)
        : levi::UnaryOperator<typename EvaluableT::row_type, EvaluableT>(expression, expression.name())
    {
        levi::unused(row);
        assert(row == 0);
        this->m_info->type = levi::EvaluableType::Row;
        this->m_info->lhs = expression;
        this->m_info->block.startRow = row;
    }

    virtual ~RowEvaluable() final;

    virtual const typename EvaluableT::row_type& evaluate() final {
        this->m_evaluationBuffer = this->m_expression.evaluate();

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename EvaluableT::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                        std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);

        return this->m_expression.getColumnDerivative(0, variable);
    }

};
template <typename EvaluableT>
levi::RowEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>::~RowEvaluable(){ }

/**
 * @brief The ColEvaluable. To retrieve a specified column from an evaluable. Specialization for matrix valued evaluables.
 */
template <typename EvaluableT>
class levi::ColEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename EvaluableT::col_type, EvaluableT>
{
    Eigen::Index m_col;

public:

    ColEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index col)
        : levi::UnaryOperator<typename EvaluableT::col_type, EvaluableT>(expression, expression.rows(), 1, "(" + expression.name() + ")(:," + std::to_string(col) + ")")
        , m_col(col)
    {
        this->m_info->type = levi::EvaluableType::Column;
        this->m_info->lhs = expression;
        this->m_info->block.startCol = col;
    }

    virtual ~ColEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<typename EvaluableT::col_type>::row_type>> row(Eigen::Index row) final {
        return this->m_expression(row, m_col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<typename EvaluableT::col_type>::value_type>> element(Eigen::Index row, Eigen::Index col) {
        levi::unused(col);
        return this->m_expression(row, m_col);
    }

    virtual const typename EvaluableT::col_type& evaluate() final {
        this->m_evaluationBuffer.lazyAssign(this->m_expression.evaluate().col(m_col));

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::col_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                            std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);
        levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::col_type>::derivative_evaluable> newDerivative;

        newDerivative = this->m_expression.getColumnDerivative(m_col, variable);

        return newDerivative;
    }

};
template <typename EvaluableT>
levi::ColEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>::~ColEvaluable(){ }

/**
 * @brief The ColEvaluable. To retrieve a specified column from an evaluable. Specialization for scalar valued evaluables.
 */
template <typename EvaluableT>
class levi::ColEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename EvaluableT::col_type, EvaluableT>
{

public:

    ColEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index col)
        : levi::UnaryOperator<typename EvaluableT::col_type, EvaluableT>(expression, expression.name())
    {
        levi::unused(col);
        assert(col == 0);
        this->m_info->type = levi::EvaluableType::Column;
        this->m_info->lhs = expression;
        this->m_info->block.startCol = col;
    }

    virtual ~ColEvaluable() final;

    virtual const typename EvaluableT::col_type& evaluate() final {
        this->m_evaluationBuffer = this->m_expression.evaluate();

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename EvaluableT::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                        std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);

        return this->m_expression.getColumnDerivative(0, variable);
    }

};
template <typename EvaluableT>
levi::ColEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>::~ColEvaluable(){ }


/**
 * @brief The ElementEvaluable. To retrieve a specified element from an evaluable. Specialization for matrix valued evaluables.
 */
template <typename EvaluableT>
class levi::ElementEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename EvaluableT::value_type, EvaluableT>
{
    Eigen::Index m_row, m_col;

public:

    ElementEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index row, Eigen::Index col)
        : levi::UnaryOperator<typename EvaluableT::value_type, EvaluableT>(expression, "[" + expression.name() + "](" + std::to_string(row) + ", " + std::to_string(col) + ")")
        , m_row(row)
        , m_col(col)
    {
        this->m_info->type = levi::EvaluableType::Element;
        this->m_info->lhs = expression;
        this->m_info->block.startRow = row;
        this->m_info->block.startCol = col;
    }

    virtual ~ElementEvaluable() final;

    virtual const typename EvaluableT::value_type& evaluate() final {
        this->m_evaluationBuffer = this->m_expression.evaluate()(m_row, m_col);

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::value_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                              std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);
        levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::value_type>::derivative_evaluable> newDerivative;

        newDerivative = levi::ExpressionComponent<levi::RowEvaluable<typename EvaluableT::derivative_evaluable>>(this->m_expression.getColumnDerivative(m_col, variable), m_row);

        return newDerivative;
    }

};
template <typename EvaluableT>
levi::ElementEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>::~ElementEvaluable(){ }

/**
 * @brief The ElementEvaluable. To retrieve a specified element from an evaluable. Specialization for scalar valued evaluables.
 */
template <typename EvaluableT>
class levi::ElementEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename EvaluableT::value_type, EvaluableT> {

public:

    ElementEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, Eigen::Index row, Eigen::Index col)
        : levi::UnaryOperator<typename EvaluableT::value_type, EvaluableT>(expression, expression.name())
    {
        levi::unused(row, col);
        assert(row == 0 && col == 0);
        this->m_info->type = levi::EvaluableType::Element;
        this->m_info->lhs = expression;
        this->m_info->block.startRow = row;
        this->m_info->block.startCol = col;
    }

    virtual ~ElementEvaluable() final;

    virtual const typename EvaluableT::value_type& evaluate() final {
        this->m_evaluationBuffer = this->m_expression.evaluate();
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename EvaluableT::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                        std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);

        return this->m_expression.getColumnDerivative(0, variable);
    }
};
template <typename EvaluableT>
levi::ElementEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>::~ElementEvaluable() { }


/**
 * @brief The BlockEvaluable. To retrieve a specified block from an evaluable. Specialization for matrix valued evaluables.
 */
template <class InputEvaluable, class OutputEvaluable>
class levi::BlockEvaluable<InputEvaluable, OutputEvaluable, typename std::enable_if<!std::is_arithmetic<typename InputEvaluable::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename OutputEvaluable::matrix_type, InputEvaluable>
{
    Eigen::Index m_startRow, m_startCol;

    template<bool value>
    void eval(levi::bool_value<value>);

    void eval(levi::bool_value<true>) {
        this->m_evaluationBuffer.lazyAssign((this->m_expression.evaluate()).template block<OutputEvaluable::rows_at_compile_time, OutputEvaluable::cols_at_compile_time>(m_startRow, m_startCol));
    }

    void eval(levi::bool_value<false>) {
        this->m_evaluationBuffer.lazyAssign(this->m_expression.evaluate().block(m_startRow, m_startCol, this->rows(), this->cols()));
    }

public:

    typedef typename OutputEvaluable::matrix_type block_type;

    BlockEvaluable(const levi::ExpressionComponent<InputEvaluable>& expression, Eigen::Index startRow, Eigen::Index startCol, Eigen::Index numberOfRows, Eigen::Index numberOfCols)
        : levi::UnaryOperator<block_type, InputEvaluable>(expression, numberOfRows, numberOfCols,
                                                          "[" + expression.name() + "](" + std::to_string(startRow) + ":" + std::to_string(startRow + numberOfRows) + ", " + std::to_string(startCol) + ":" + std::to_string(startCol + numberOfCols) + ")")
        , m_startRow(startRow)
        , m_startCol(startCol)
    {
        assert(((startRow + numberOfRows) <= expression.rows()) && ((startCol + numberOfCols) <= expression.cols()));
        this->m_info->type = levi::EvaluableType::Block;
        this->m_info->lhs = expression;
        this->m_info->block.startRow = startRow;
        this->m_info->block.startCol = startCol;
        this->m_info->block.rows = numberOfRows;
        this->m_info->block.cols = numberOfCols;
    }

    virtual ~BlockEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<block_type>::row_type>> row(Eigen::Index row) final {
        return this->m_expression.row(m_startRow + row).block(0, m_startCol, 1, this->cols());
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<block_type>::col_type>> col(Eigen::Index col) final {
        return this->m_expression.col(m_startCol + col).block(m_startRow, 0, this->rows(), 1);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<block_type>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return this->m_expression(m_startRow + row, m_startCol + col);
    }

    virtual const block_type& evaluate() final {

        eval(levi::bool_value<OutputEvaluable::rows_at_compile_time != Eigen::Dynamic && OutputEvaluable::cols_at_compile_time != Eigen::Dynamic>());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<block_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {
        levi::ExpressionComponent<typename levi::Evaluable<block_type>::derivative_evaluable> newDerivative;

        newDerivative = this->m_expression.getColumnDerivative(m_startCol + column, variable).block(m_startRow, 0, this->rows(), variable->dimension());

        return newDerivative;
    }

};
template <class InputEvaluable, class OutputEvaluable>
levi::BlockEvaluable<InputEvaluable, OutputEvaluable, typename std::enable_if<!std::is_arithmetic<typename InputEvaluable::matrix_type>::value>::type>::~BlockEvaluable(){ }

/**
 * @brief The BlockEvaluable. To retrieve a specified block from an evaluable. Specialization for scalar valued evaluables.
 */
template <class InputEvaluable, class OutputEvaluable>
class levi::BlockEvaluable<InputEvaluable, OutputEvaluable, typename std::enable_if<std::is_arithmetic<typename InputEvaluable::matrix_type>::value>::type>
        : public levi::UnaryOperator<typename OutputEvaluable::matrix_type, InputEvaluable> {

public:

    typedef typename OutputEvaluable::matrix_type block_type;

    BlockEvaluable(const levi::ExpressionComponent<InputEvaluable>& expression, Eigen::Index startRow, Eigen::Index startCol, Eigen::Index numberOfRows, Eigen::Index numberOfCols)
        : levi::UnaryOperator<block_type, InputEvaluable>(expression, expression.name())
    {
        levi::unused(startRow, startCol, numberOfRows, numberOfCols);
        assert(startRow == 0 && startCol == 0 && numberOfRows == 1 && numberOfCols == 1);
        this->m_info->type = levi::EvaluableType::Block;
        this->m_info->lhs = expression;
        this->m_info->block.startRow = startRow;
        this->m_info->block.startCol = startCol;
        this->m_info->block.rows = numberOfRows;
        this->m_info->block.cols = numberOfCols;
    }

    virtual ~BlockEvaluable() final;

    virtual const block_type& evaluate() final {
        this->m_evaluationBuffer = this->m_expression.evaluate();
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<block_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                         std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        assert(column == 0);

        return this->m_expression.getColumnDerivative(0, variable);
    }

};
template <class InputEvaluable, class OutputEvaluable>
levi::BlockEvaluable<InputEvaluable, OutputEvaluable, typename std::enable_if<std::is_arithmetic<typename InputEvaluable::matrix_type>::value>::type>::~BlockEvaluable(){ }


/**
 * @brief The CastEvaluable.
 *
 * It allows to assign an evaluable to an expression whose pointer cannot be directly casted. It assumes that the two evaluation buffers can be casted.
 */
template <class LeftEvaluable, class RightEvaluable>
class levi::CastEvaluable : public levi::UnaryOperator<typename LeftEvaluable::matrix_type, RightEvaluable> {

    template <bool leftIsScalar, bool rightIsScalar>
    void eval(levi::bool_value<leftIsScalar>, levi::bool_value<rightIsScalar>());

    void eval(levi::bool_value<true>, levi::bool_value<true>) {
        this->m_evaluationBuffer = this->m_expression.evaluate();
    }

    void eval(levi::bool_value<true>, levi::bool_value<false>) {
        this->m_evaluationBuffer = this->m_expression.evaluate()(0,0);
    }

    void eval(levi::bool_value<false>, levi::bool_value<true>) {
        this->m_evaluationBuffer(0,0) = this->m_expression.evaluate();
    }

    template <bool ImpossibleAssignement>
    void evalMatrix(levi::bool_value<ImpossibleAssignement>());

    void evalMatrix(levi::bool_value<true>) {
        assert(false && "You tried to assign Evaluables of different sizes.");
    }

    void evalMatrix(levi::bool_value<false>) {
        this->m_evaluationBuffer.lazyAssign(this->m_expression.evaluate());
    }

    void eval(levi::bool_value<false>, levi::bool_value<false>) {
        evalMatrix(levi::bool_value<((LeftEvaluable::rows_at_compile_time * RightEvaluable::rows_at_compile_time > 0) && (LeftEvaluable::rows_at_compile_time != RightEvaluable::rows_at_compile_time)) ||
                   ((LeftEvaluable::cols_at_compile_time * RightEvaluable::cols_at_compile_time > 0) && (LeftEvaluable::cols_at_compile_time != RightEvaluable::cols_at_compile_time))>());
    }

public:

    CastEvaluable(const levi::ExpressionComponent<RightEvaluable>& rhs)
        : levi::UnaryOperator<typename LeftEvaluable::matrix_type, RightEvaluable>(rhs, rhs.rows(), rhs.cols(), rhs.name())
    {
        this->m_info->copy(rhs.info());
    }

    virtual ~CastEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename LeftEvaluable::row_type>> row(Eigen::Index row) final {
        return this->m_expression.row(row);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename LeftEvaluable::col_type>> col(Eigen::Index col) final {
        return this->m_expression.col(col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename LeftEvaluable::value_type>> element(Eigen::Index row, Eigen::Index col) {
        return this->m_expression(row, col);
    }

    virtual const typename LeftEvaluable::matrix_type& evaluate() final {
        eval(levi::bool_value<std::is_arithmetic<typename LeftEvaluable::matrix_type>::value>(), levi::bool_value<std::is_arithmetic<typename RightEvaluable::matrix_type>::value>());
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename LeftEvaluable::matrix_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::ExpressionComponent<typename RightEvaluable::derivative_evaluable> rightDerivative = this->m_expression.getColumnDerivative(column, variable);

        levi::ExpressionComponent<CastEvaluable<typename levi::Evaluable<typename LeftEvaluable::matrix_type>::derivative_evaluable, typename RightEvaluable::derivative_evaluable>> newCast(rightDerivative);

        return newCast;
    }
};
template <class LeftEvaluable, class RightEvaluable>
levi::CastEvaluable<LeftEvaluable, RightEvaluable>::~CastEvaluable() {}

#endif // LEVI_ACCESSOREVALUABLES_H
