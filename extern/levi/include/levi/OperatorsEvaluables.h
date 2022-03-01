/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_OPERATORS_EVALUABLES_H
#define LEVI_OPERATORS_EVALUABLES_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/BasicEvaluables.h>
#include <levi/Expression.h>
#include <levi/OperatorsBase.h>

/**
 * Helper struct for determining the type resulting from an addition
 */
template <typename Scalar_lhs, typename Scalar_rhs>
struct levi::scalar_sum_return {
    //decltype allow to get the return type of the addition of a variable of type Scalar_lhs to a variable of type Scalar_rhs.
    typedef decltype (std::declval<Scalar_lhs>() + std::declval<Scalar_rhs>()) type;
};

/**
 * Helper struct for determining the type resulting from an addition of two matrices. Specialization for two matrices
 *
 */
template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_sum_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
        typename std::enable_if<levi::is_valid_sum<lhsRows, lhsCols, rhsRows, rhsCols>::value>::type> {
    typedef Eigen::Matrix<typename levi::scalar_sum_return<Scalar_lhs, Scalar_rhs>::type, (lhsRows < rhsRows) ? lhsRows : rhsRows, (lhsCols < rhsCols) ? lhsCols : rhsCols> type; //here we assume that Eigen::Dynamic == -1, thus, given that it is a valid sum, the minimum will be -1 if present, thus using Eigen::Dynamic as output
};

/**
 * Helper struct for determining the type resulting from an addition of two matrices. Specialization for two scalars.
 *
 */
template<typename Scalar_lhs, typename Scalar_rhs>
struct levi::matrix_sum_return<Scalar_lhs, Scalar_rhs,
        typename std::enable_if<std::is_arithmetic<Scalar_lhs>::value && std::is_arithmetic<Scalar_rhs>::value>::type> {
    typedef typename levi::scalar_sum_return<Scalar_lhs, Scalar_rhs>::type type;
};

/**
 * Helper struct for determining the type resulting from an addition of two matrices. Specialization for a scalar and a matrix.
 *
 */
template<typename Scalar, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_sum_return<Scalar, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
        typename std::enable_if<std::is_arithmetic<Scalar>::value && levi::is_valid_sum<1,1, rhsRows, rhsCols>::value>::type> {
    typedef Eigen::Matrix<typename levi::scalar_sum_return<Scalar, Scalar_rhs>::type, rhsRows, rhsCols> type;
};

/**
 * Helper struct for determining the type resulting from an addition of two matrices. Specialization for a matrix and a scalar
 *
 */
template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar>
struct levi::matrix_sum_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Scalar,
        typename std::enable_if<std::is_arithmetic<Scalar>::value && levi::is_valid_sum<1,1, lhsRows, lhsCols>::value>::type> {
    typedef Eigen::Matrix<typename levi::scalar_sum_return<Scalar, Scalar_lhs>::type, lhsRows, lhsCols> type;
};

/**
 * Helper struct for determining the type resulting from a multiplication
 */
template <typename Scalar_lhs, typename Scalar_rhs>
struct levi::scalar_product_return {
    //decltype allow to get the return type of the multiplication of a variable of type Scalar_lhs to a variable of type Scalar_rhs.
    typedef decltype (std::declval<Scalar_lhs>() * std::declval<Scalar_rhs>()) type;
};

/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization for two matrices
 *
 */
template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
        typename std::enable_if<levi::is_valid_product<lhsRows, lhsCols, rhsRows, rhsCols>::value && !(lhsRows == 1 && lhsCols == 1) && !(rhsRows == 1 && rhsCols == 1)>::type> {
    typedef Eigen::Matrix<typename levi::scalar_product_return<Scalar_lhs, Scalar_rhs>::type, lhsRows, rhsCols> type;
};

/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization in the case the lhs has exactly one row and one column.
 *
 */
template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
        typename std::enable_if<lhsRows == 1 && lhsCols == 1 && (rhsRows != 1 || rhsCols != 1)>::type> {
    typedef Eigen::Matrix<typename levi::scalar_product_return<Scalar_lhs, Scalar_rhs>::type, rhsRows, rhsCols> type;
};


/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization in the case the rhs has exactly one row and one column.
 *
 */
template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
        typename std::enable_if<rhsRows == 1 && rhsCols == 1 && (lhsRows != 1 || lhsCols != 1)>::type> {
    typedef Eigen::Matrix<typename levi::scalar_product_return<Scalar_lhs, Scalar_rhs>::type, lhsRows, lhsCols> type;
};

/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization in the case both the rhs and the lhs have exactly one row and one column.
 *
 */
template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
                                   typename std::enable_if<rhsRows == 1 && rhsCols == 1 && lhsRows == 1 && lhsCols == 1>::type> {
    typedef Eigen::Matrix<typename levi::scalar_product_return<Scalar_lhs, Scalar_rhs>::type, 1, 1> type;
};

/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization for a scalar and a matrix
 */
template<typename Scalar, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
struct levi::matrix_product_return<Scalar, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
        typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> {
    typedef Eigen::Matrix<typename levi::scalar_product_return<Scalar, Scalar_rhs>::type, rhsRows, rhsCols> type;
};

/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization for a matrix and a scalar.
 */
template<typename Scalar, typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols>
struct levi::matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Scalar,
        typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> {
    typedef Eigen::Matrix<typename levi::scalar_product_return<Scalar, Scalar_lhs>::type, lhsRows, lhsCols> type;
};

/**
 * Helper struct for determining the type resulting from a multiplication of two matrices. Specialization for two scalars.
 */
template<typename Scalar_lhs, typename Scalar_rhs>
struct levi::matrix_product_return<Scalar_lhs, Scalar_rhs,
        typename std::enable_if<std::is_arithmetic<Scalar_lhs>::value && std::is_arithmetic<Scalar_rhs>::value>::type> {
    typedef typename levi::scalar_product_return<Scalar_lhs, Scalar_rhs>::type type;
};

/**
 * @brief The SumEvaluable. Implements the sum of two evaluables.
 */
template <class LeftEvaluable, class RightEvaluable>
class levi::SumEvaluable : public levi::BinaryOperator<typename levi::matrix_sum_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::matrix_type>::type, LeftEvaluable, RightEvaluable>{

    template<bool value>
    void eval(levi::bool_value<value>);

    void eval(levi::bool_value<true>) {
        this->m_evaluationBuffer = this->m_lhs.evaluate() + this->m_rhs.evaluate();
    }

    void eval(levi::bool_value<false>) {
        this->m_evaluationBuffer.lazyAssign(this->m_lhs.evaluate() + this->m_rhs.evaluate());
    }

public:

    typedef typename levi::matrix_sum_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::matrix_type>::type sum_type;

    SumEvaluable(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs)
        : levi::BinaryOperator<sum_type, LeftEvaluable, RightEvaluable>(lhs, rhs, lhs.rows(), lhs.cols(), "(" + lhs.name() + " + " + rhs.name() + ")")
    {
        this->m_info->type = levi::EvaluableType::Sum;
        this->m_info->lhs = lhs;
        this->m_info->rhs = rhs;
    }

    virtual ~SumEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<sum_type>::row_type>> row(Eigen::Index row) final {
        return this->m_lhs.row(row) + this->m_rhs.row(row);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<sum_type>::col_type>> col(Eigen::Index col) final {
        return this->m_lhs.col(col) + this->m_rhs.col(col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<sum_type>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return this->m_lhs(row, col) + this->m_rhs(row, col);
    }

    virtual const sum_type& evaluate() final {

        eval(levi::bool_value<std::is_arithmetic<sum_type>::value>());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<sum_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                       std::shared_ptr<levi::VariableBase> variable) final {
        levi::ExpressionComponent<typename levi::Evaluable<sum_type>::derivative_evaluable> sumDerivative;

        bool isLeftDependent = this->m_lhs.isDependentFrom(variable);
        bool isRightDependent = this->m_rhs.isDependentFrom(variable);

        if (isLeftDependent && isRightDependent) {
            sumDerivative = this->m_lhs.getColumnDerivative(column, variable) + this->m_rhs.getColumnDerivative(column, variable);

            return sumDerivative;
        }

        if (!isLeftDependent && !isRightDependent) {
            levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<sum_type>::derivative_evaluable::matrix_type>> nullDerivative(this->rows(), variable->dimension(), "d" + this->name() + "/d" + variable->variableName());
            sumDerivative = nullDerivative;
            return sumDerivative;
        }

        if (isLeftDependent) {
            sumDerivative = this->m_lhs.getColumnDerivative(column, variable);
            return sumDerivative;
        } else {
            sumDerivative = this->m_rhs.getColumnDerivative(column, variable);
            return sumDerivative;
        }
    }

};
template <class LeftEvaluable, class RightEvaluable>
levi::SumEvaluable<LeftEvaluable, RightEvaluable>::~SumEvaluable() { }

/**
 * @brief The SubtractionEvaluable. Implements the subtraction of two evaluables.
 */
template <class LeftEvaluable, class RightEvaluable>
class levi::SubtractionEvaluable : public levi::BinaryOperator<typename levi::matrix_sum_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::matrix_type>::type, LeftEvaluable, RightEvaluable>{

    template<bool value>
    void eval(levi::bool_value<value>);

    void eval(levi::bool_value<true>) {
        this->m_evaluationBuffer = this->m_lhs.evaluate() - this->m_rhs.evaluate();
    }

    void eval(levi::bool_value<false>) {
        this->m_evaluationBuffer.lazyAssign(this->m_lhs.evaluate() - this->m_rhs.evaluate());
    }

public:

    typedef typename levi::matrix_sum_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::matrix_type>::type sum_type;

    SubtractionEvaluable(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs)
        : levi::BinaryOperator<sum_type, LeftEvaluable, RightEvaluable>(lhs, rhs, lhs.rows(), lhs.cols(), "(" + lhs.name() + " - " + rhs.name() + ")")
    {
        this->m_info->type = levi::EvaluableType::Subtraction;
        this->m_info->lhs = lhs;
        this->m_info->rhs = rhs;
    }

    virtual ~SubtractionEvaluable() final;

    virtual const sum_type& evaluate() final {
        eval(levi::bool_value<std::is_arithmetic<sum_type>::value>());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<sum_type>::row_type>> row(Eigen::Index row) final {
        return this->m_lhs.row(row) - this->m_rhs.row(row);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<sum_type>::col_type>> col(Eigen::Index col) final {
        return this->m_lhs.col(col) - this->m_rhs.col(col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<sum_type>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return this->m_lhs(row, col) - this->m_rhs(row, col);
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<sum_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                       std::shared_ptr<levi::VariableBase> variable) final {
        levi::ExpressionComponent<typename levi::Evaluable<sum_type>::derivative_evaluable> sumDerivative;

        bool isLeftDependent = this->m_lhs.isDependentFrom(variable);
        bool isRightDependent = this->m_rhs.isDependentFrom(variable);

        if (isLeftDependent && isRightDependent) {
            sumDerivative = this->m_lhs.getColumnDerivative(column, variable) - this->m_rhs.getColumnDerivative(column, variable);

            return sumDerivative;
        }

        if (!isLeftDependent && !isRightDependent) {
            levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<sum_type>::derivative_evaluable::matrix_type>> nullDerivative(this->rows(), variable->dimension(), "d" + this->name() + "/d" + variable->variableName());
            sumDerivative = nullDerivative;
            return sumDerivative;
        }

        if (isLeftDependent) {
            sumDerivative = this->m_lhs.getColumnDerivative(column, variable);
            return sumDerivative;
        } else {
            sumDerivative = - this->m_rhs.getColumnDerivative(column, variable);
            return sumDerivative;
        }
    }
};
template <class LeftEvaluable, class RightEvaluable>
levi::SubtractionEvaluable<LeftEvaluable, RightEvaluable>::~SubtractionEvaluable() { }

template <class EvaluableT>
class levi::SignInvertedEvaluable : public levi::UnaryOperator<typename EvaluableT::matrix_type, EvaluableT>{

    template<bool value>
    void eval(levi::bool_value<value>);

    void eval(levi::bool_value<true>) {
        this->m_evaluationBuffer = -(this->m_expression.evaluate());
    }

    void eval(levi::bool_value<false>) {
        this->m_evaluationBuffer.lazyAssign(-(this->m_expression.evaluate()));
    }

public:

    SignInvertedEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, int)
        : levi::UnaryOperator<typename EvaluableT::matrix_type, EvaluableT>(expression, expression.rows(), expression.cols(), "-" + expression.name())
    {
        this->m_info->lhs = expression;
        this->m_info->type = levi::EvaluableType::InvertedSign;
    }

    virtual ~SignInvertedEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename EvaluableT::row_type>> row(Eigen::Index row) final {
        return -this->m_expression.row(row);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename EvaluableT::col_type>> col(Eigen::Index col) final {
        return -this->m_expression.col(col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename EvaluableT::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return -this->m_expression(row, col);
    }

    virtual const typename EvaluableT::matrix_type& evaluate() final {
        eval(levi::bool_value<std::is_arithmetic<typename EvaluableT::matrix_type>::value>());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::matrix_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                               std::shared_ptr<levi::VariableBase> variable) final {
        levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::matrix_type>::derivative_evaluable> derivative;

        derivative = -(this->m_expression.getColumnDerivative(column, variable));

        return derivative;
    }
};
template <class EvaluableT>
levi::SignInvertedEvaluable<EvaluableT>::~SignInvertedEvaluable() { }

template <class LeftEvaluable, class RightEvaluable>
class levi::MatrixProductDerivative : public levi::Evaluable<
                                          Eigen::Matrix<typename levi::scalar_product_return<typename LeftEvaluable::value_type,
                                                                                             typename RightEvaluable::value_type>::type,
                                                        LeftEvaluable::rows_at_compile_time, Eigen::Dynamic>> {

    levi::ExpressionComponent<LeftEvaluable> m_lhs;
    levi::ExpressionComponent<RightEvaluable> m_rhs;

    using return_type = Eigen::Matrix<typename levi::scalar_product_return<typename LeftEvaluable::value_type,
                                                                           typename RightEvaluable::value_type>::type,
                                      LeftEvaluable::rows_at_compile_time, Eigen::Dynamic>;

    using this_derivative_evaluable = typename levi::Evaluable<return_type>::derivative_evaluable;

    using lhs_derivative_evaluable = typename LeftEvaluable::derivative_evaluable;

    using lhs_derivative_evaluable_cols = typename lhs_derivative_evaluable::col_type;

    using composite_matrix = Eigen::Matrix<typename LeftEvaluable::value_type, LeftEvaluable::rows_at_compile_time, Eigen::Dynamic>;

    using single_variable_derivative_expression =
        levi::ExpressionComponent<levi::Evaluable<composite_matrix>>;

    std::vector<levi::ExpressionComponent<lhs_derivative_evaluable>> m_lhsDerivatives;

    std::vector<single_variable_derivative_expression> m_lhsSingleVariableDerivatives;

public:


    MatrixProductDerivative(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs, std::shared_ptr<levi::VariableBase> variable)
        : levi::Evaluable<return_type>(lhs.rows(), variable->dimension(),
                                       "d(" + lhs.name() + ")/d" + variable->variableName() + " * " + rhs.name())
          , m_lhs(lhs)
          , m_rhs(rhs)
    {
        static_assert (!std::is_arithmetic<typename LeftEvaluable::matrix_type>::value, "The left hand side is supposed to be a matrix");
        static_assert (RightEvaluable::cols_at_compile_time == 1 || RightEvaluable::cols_at_compile_time == Eigen::Dynamic,
                      "The right hand side is supposed to be a vector");

        assert(lhs.isValidExpression() && "The left hand side does not appear to be a valid expression.");
        assert(rhs.isValidExpression() && "The right hand side does not appear to be a valid expression.");
        assert(lhs.cols() == rhs.rows() && "The left hand side columns do not match the right hand side rows.");

        this->m_info->lhs = lhs;
        this->m_info->rhs = rhs;

        this->addDependencies(m_rhs);

        for (Eigen::Index col = 0; col < lhs.cols(); ++col) {
            m_lhsDerivatives.push_back(m_lhs.getColumnDerivative(col, variable));
            this->addDependencies(m_lhsDerivatives.back());
        }

        for (Eigen::Index variableDim = 0; variableDim < variable->dimension(); ++variableDim) {
            std::vector<levi::ExpressionComponent<levi::Evaluable<lhs_derivative_evaluable_cols>>> derivative_cols;
            for (Eigen::Index col = 0; col < lhs.cols(); ++col) {
                derivative_cols.push_back(m_lhsDerivatives[col].col(variableDim));
            }
            m_lhsSingleVariableDerivatives.push_back(levi::ComposeByCols(derivative_cols, "d(" + lhs.name() + ")/d" +
                                                                             variable->variableName() + "_" +
                                                                             std::to_string(variableDim)));
        }

    }

    virtual ~MatrixProductDerivative() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<return_type>::col_type>> col(Eigen::Index col) final {
        return m_lhsSingleVariableDerivatives[col] * m_rhs;
    }

    virtual const return_type& evaluate() final {

        const typename RightEvaluable::matrix_type& rhs = m_rhs.evaluate();

        this->m_evaluationBuffer =rhs(0,0) * m_lhsDerivatives[0].evaluate();

        for (Eigen::Index col = 1; col < m_lhsDerivatives.size(); ++col) {
            this->m_evaluationBuffer += rhs(col,0) * m_lhsDerivatives[col].evaluate();
        }

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<this_derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column,
                           std::shared_ptr<levi::VariableBase> variable) final {

                return col(column).getColumnDerivative(0, variable);
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_lhs.clearDerivativesCache();
        m_rhs.clearDerivativesCache();

        for (auto& expression : m_lhsDerivatives) {
            expression.clearDerivativesCache();
        }

        for (auto& expression : m_lhsSingleVariableDerivatives) {
            expression.clearDerivativesCache();
        }
    }
};
template <class LeftEvaluable, class RightEvaluable>
levi::MatrixProductDerivative<LeftEvaluable, RightEvaluable>::~MatrixProductDerivative() { }

/**
 * @brief The ProductEvaluable. Implements the product of two evaluables.
 */
template <class LeftEvaluable, class RightEvaluable>
class levi::ProductEvaluable : public levi::BinaryOperator<typename levi::matrix_product_return<
        typename LeftEvaluable::matrix_type, typename RightEvaluable::matrix_type>::type, LeftEvaluable, RightEvaluable>{

public:

    typedef typename levi::matrix_product_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::matrix_type>::type product_type;

private:

    typedef levi::ExpressionComponent<typename levi::Evaluable<product_type>::derivative_evaluable> derivative_expression;

    template<typename Matrix>
    inline const Matrix& cast_unit_matrix(levi::bool_value<false>, const Matrix& input) {
        return input;
    }

    template<typename Matrix>
    inline Matrix& cast_unit_matrix(levi::bool_value<false>, Matrix& input) {
        return input;
    }

    template<typename Scalar>
    inline const Scalar& cast_unit_matrix(levi::bool_value<true>, const Eigen::Matrix<Scalar, 1,1>& input) {
        return input(0,0);
    }

    template<typename Scalar>
    inline Scalar& cast_unit_matrix(levi::bool_value<true>, Eigen::Matrix<Scalar, 1,1>& input) {
        return input(0,0);
    }

    template<bool lhsIsScalar, bool rhsIsScalar>
    derivative_expression get_derivative(levi::bool_value<lhsIsScalar>, levi::bool_value<rhsIsScalar>, Eigen::Index column, std::shared_ptr<levi::VariableBase> variable);

    /**
     * @brief Helper function for the derivative of the multiplication between two scalars
     */
    derivative_expression get_derivative(levi::bool_value<true>, levi::bool_value<true>, Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
        derivative_expression derivative;

        bool isLeftDependent = this->m_lhs.isDependentFrom(variable);
        bool isRightDependent = this->m_rhs.isDependentFrom(variable);

        if (isLeftDependent && isRightDependent) {
            derivative = this->m_rhs * this->m_lhs.getColumnDerivative(column, variable) + this->m_lhs * this->m_rhs.getColumnDerivative(column, variable);
            return derivative;
        }

        if (!isLeftDependent && !isRightDependent) {
            levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<product_type>::derivative_evaluable::matrix_type>> nullDerivative(this->rows(), variable->dimension(), "d(" + this->name() + ")/d" + variable->variableName());
            derivative = nullDerivative;
            return derivative;
        }

        if (isLeftDependent) {
            derivative = this->m_rhs * this->m_lhs.getColumnDerivative(column, variable);
            return derivative;
        } else {
            derivative = this->m_lhs * this->m_rhs.getColumnDerivative(column, variable);
            return derivative;
        }
    }

    /**
     * @brief Helper function for the derivative of the multiplication between a scalar and a matrix.
     */
    derivative_expression get_derivative(levi::bool_value<true>, levi::bool_value<false>, Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {

        //k*A
        derivative_expression derivative;

        bool isLeftDependent = this->m_lhs.isDependentFrom(variable);
        bool isRightDependent = this->m_rhs.isDependentFrom(variable);

        if (isLeftDependent && isRightDependent) {
            derivative = this->m_lhs * this->m_rhs.getColumnDerivative(column, variable) + this->m_rhs.col(column) * this->m_lhs.getColumnDerivative(0, variable);
            return derivative;
        }

        if (!isLeftDependent && !isRightDependent) {
            levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<product_type>::derivative_evaluable::matrix_type>> nullDerivative(this->rows(), variable->dimension(), "d(" + this->name() + ")/d" + variable->variableName());
            derivative = nullDerivative;
            return derivative;
        }

        if (isLeftDependent) {
            derivative = this->m_rhs.col(column) * this->m_lhs.getColumnDerivative(0, variable);
            return derivative;
        } else {
            derivative = this->m_lhs * this->m_rhs.getColumnDerivative(column, variable);
            return derivative;
        }
    }

    /**
     * @brief Helper function for the derivative of the multiplication between a matrix and a scalar.
     */
    derivative_expression get_derivative(levi::bool_value<false>, levi::bool_value<true>, Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {

        //A*k
        derivative_expression derivative;

        bool isLeftDependent = this->m_lhs.isDependentFrom(variable);
        bool isRightDependent = this->m_rhs.isDependentFrom(variable);

        if (isLeftDependent && isRightDependent) {
            derivative = this->m_rhs * this->m_lhs.getColumnDerivative(column, variable) + this->m_lhs.col(column) * this->m_rhs.getColumnDerivative(0, variable);
            return derivative;
        }

        if (!isLeftDependent && !isRightDependent) {
            levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<product_type>::derivative_evaluable::matrix_type>> nullDerivative(this->rows(), variable->dimension(), "d(" + this->name() + ")/d" + variable->variableName());
            derivative = nullDerivative;
            return derivative;
        }

        if (isLeftDependent) {
            derivative = this->m_rhs * this->m_lhs.getColumnDerivative(column, variable);
            return derivative;
        } else {
            derivative = this->m_lhs.col(column) * this->m_rhs.getColumnDerivative(0, variable);
            return derivative;
        }
    }

    /**
     * @brief Helper function for the derivative of the multiplication between two matrices.
     */
    derivative_expression get_derivative(levi::bool_value<false>, levi::bool_value<false>, Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
        derivative_expression derivative;

        bool isLeftDependent = this->m_lhs.isDependentFrom(variable);
        bool isRightDependent = this->m_rhs.isDependentFrom(variable);

        if (isLeftDependent && isRightDependent) {

            derivative = this->m_lhs * this->m_rhs.getColumnDerivative(column, variable) +
                levi::MatrixProductDerivativeExpression(this->m_lhs, this->m_rhs.col(column), variable);

            return derivative;

        }

        if (!isLeftDependent && !isRightDependent) {
            levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<product_type>::derivative_evaluable::matrix_type>> nullDerivative(this->rows(), variable->dimension(), "d(" + this->name() + ")/d" + variable->variableName());
            derivative = nullDerivative;
            return derivative;
        }

        if (isLeftDependent) {

            derivative = levi::MatrixProductDerivativeExpression(this->m_lhs, this->m_rhs.col(column), variable);
            return derivative;

        } else {

            derivative = this->m_lhs * this->m_rhs.getColumnDerivative(column, variable);

            return derivative;
        }

    }

    inline void eval(levi::bool_value<true>) {
        cast_unit_matrix(levi::bool_value<!std::is_arithmetic<product_type>::value>(), this->m_evaluationBuffer) = cast_unit_matrix(levi::bool_value<!std::is_arithmetic<typename LeftEvaluable::matrix_type>::value>(), this->m_lhs.evaluate()) *
            cast_unit_matrix(levi::bool_value<!std::is_arithmetic<typename RightEvaluable::matrix_type>::value>(), this->m_rhs.evaluate());
    }

    inline void eval(levi::bool_value<false>) {
        this->m_evaluationBuffer.lazyAssign(cast_unit_matrix(levi::bool_value<!std::is_arithmetic<typename LeftEvaluable::matrix_type>::value && LeftEvaluable::rows_at_compile_time == 1 && LeftEvaluable::cols_at_compile_time == 1>(),this->m_lhs.evaluate()) *
                                            cast_unit_matrix(levi::bool_value<!std::is_arithmetic<typename RightEvaluable::matrix_type>::value && RightEvaluable::rows_at_compile_time == 1 && RightEvaluable::cols_at_compile_time == 1>(),this->m_rhs.evaluate()));
    }

    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::row_type>> getRow(levi::bool_value<true>, Eigen::Index row) {
        return this->m_lhs * this->m_rhs.row(row);
    }

    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::row_type>> getRow(levi::bool_value<false>, Eigen::Index row) {
        return this->m_lhs.row(row) * this->m_rhs;
    }

    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::col_type>> getCol(levi::bool_value<true>, Eigen::Index col) {
        return this->m_lhs.col(col) * this->m_rhs;
    }

    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::col_type>> getCol(levi::bool_value<false>, Eigen::Index col) {
        return this->m_lhs * this->m_rhs.col(col);
    }

    //Both operands are matrices
    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::value_type>> getElement(levi::bool_value<false>,
                                                                                                                     levi::bool_value<false>,
                                                                                                                     Eigen::Index row,
                                                                                                                     Eigen::Index col) {
        return this->m_lhs.row(row) * this->m_rhs.col(col);
    }

    //Left operand is either a scalar or a unit matrix
    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::value_type>> getElement(levi::bool_value<true>,
                                                                                                                     levi::bool_value<false>,
                                                                                                                     Eigen::Index row,
                                                                                                                     Eigen::Index col) {
        return this->m_lhs * this->m_rhs(row, col);
    }

    //Left operand is either a scalar or a unit matrix
    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::value_type>> getElement(levi::bool_value<false>,
                                                                                                                     levi::bool_value<true>,
                                                                                                                     Eigen::Index row,
                                                                                                                     Eigen::Index col) {
        return this->m_lhs(row, col) * this->m_rhs;
    }

    //Left operand is either a scalar or a unit matrix
    inline levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::value_type>> getElement(levi::bool_value<true>,
                                                                                                                     levi::bool_value<true>,
                                                                                                                     Eigen::Index row,
                                                                                                                     Eigen::Index col) {
        levi::unused(row, col);
        return this->m_lhs * this->m_rhs;
    }



public:


    ProductEvaluable(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs)
        : levi::BinaryOperator<product_type, LeftEvaluable, RightEvaluable>(lhs, rhs, (lhs.rows() == 1 && lhs.cols() == 1 && rhs.rows() != 1)? rhs.rows() : lhs.rows(),
                                        (rhs.rows() == 1 && rhs.cols() == 1 && lhs.cols() != 1)? lhs.cols() : rhs.cols(),
                                        lhs.name() + " * " + rhs.name())
    {
        this->m_info->type = levi::EvaluableType::Product;
        this->m_info->lhs = lhs;
        this->m_info->rhs = rhs;
    }

    virtual ~ProductEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::row_type>> row(Eigen::Index row) final {
        return getRow(levi::bool_value<std::is_arithmetic<typename LeftEvaluable::matrix_type>::value || (LeftEvaluable::rows_at_compile_time == 1 && LeftEvaluable::cols_at_compile_time == 1)>(), row);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::col_type>> col(Eigen::Index col) final {
        return getCol(levi::bool_value<std::is_arithmetic<typename RightEvaluable::matrix_type>::value || (RightEvaluable::rows_at_compile_time == 1 && RightEvaluable::cols_at_compile_time == 1)>(), col);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return getElement(levi::bool_value<std::is_arithmetic<typename LeftEvaluable::matrix_type>::value || (LeftEvaluable::rows_at_compile_time == 1 && LeftEvaluable::cols_at_compile_time == 1)>(),
                          levi::bool_value<std::is_arithmetic<typename RightEvaluable::matrix_type>::value || (RightEvaluable::rows_at_compile_time == 1 && RightEvaluable::cols_at_compile_time == 1)>(),
                          row, col);
    }

    virtual const product_type& evaluate() final {

        eval(levi::bool_value<std::is_arithmetic<product_type>::value ||
                              (LeftEvaluable::rows_at_compile_time == 1 && LeftEvaluable::cols_at_compile_time == 1
              && RightEvaluable::rows_at_compile_time == 1 && RightEvaluable::cols_at_compile_time == 1)>());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<product_type>::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column,
                           std::shared_ptr<levi::VariableBase> variable) final {

        return get_derivative(levi::bool_value<std::is_arithmetic<typename LeftEvaluable::matrix_type>::value>(),
                              levi::bool_value<std::is_arithmetic<typename RightEvaluable::matrix_type>::value>(),
                              column, variable);
    }
};
template <class LeftEvaluable, class RightEvaluable>
levi::ProductEvaluable<LeftEvaluable, RightEvaluable>::~ProductEvaluable() { }


template <class EvaluableT>
class levi::PowEvaluable : public levi::UnaryOperator<typename EvaluableT::value_type, EvaluableT> {

    typename EvaluableT::value_type m_exponent;

    template<bool rhsIsScalar>
    double get_value(levi::bool_value<rhsIsScalar>);

    double get_value(levi::bool_value<true>) {
        return this->m_expression.evaluate();
    }

    double get_value(levi::bool_value<false>) {
        return this->m_expression.evaluate()(0,0);
    }

public:

    PowEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, typename EvaluableT::value_type exponent)
        : levi::UnaryOperator<typename EvaluableT::value_type, EvaluableT>(expression, "(" +expression.name() + ")^(" + std::to_string(exponent) + ")")
        , m_exponent(exponent)
    {
        assert(this->m_expression.rows() == 1 && this->m_expression.cols() == 1);
        this->m_info->type = levi::EvaluableType::Pow;
        this->m_info->exponent = exponent;
        this->m_info->lhs = expression;
    }

    virtual ~PowEvaluable() final;

    virtual const typename EvaluableT::value_type& evaluate() final {

        this->m_evaluationBuffer = std::pow(get_value(levi::bool_value<std::is_arithmetic<typename EvaluableT::matrix_type>::value>()), m_exponent);

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::value_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                              std::shared_ptr<levi::VariableBase> variable) final {
        assert(column == 0);

        levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::value_type>::derivative_evaluable> derivative;

        if (std::abs(m_exponent -1) < 1e-15) {
            derivative = this->m_expression.getColumnDerivative(column, variable);
        } else {
            derivative = this->m_expression.pow(m_exponent - 1) * m_exponent * this->m_expression.getColumnDerivative(column, variable);
        }

        return derivative;
    }
};
template <class EvaluableT>
levi::PowEvaluable<EvaluableT>::~PowEvaluable() { }

template <class LeftEvaluable, class RightEvaluable>
class levi::DivisionEvaluable : public levi::BinaryOperator<typename levi::matrix_product_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::value_type>::type, LeftEvaluable, RightEvaluable>
{
public:

    typedef typename levi::matrix_product_return<typename LeftEvaluable::matrix_type, typename RightEvaluable::value_type>::type product_type;

private:

    template<bool rhsIsScalar>
    double get_rhs_value(levi::bool_value<rhsIsScalar>);

    double get_rhs_value(levi::bool_value<true>) {
        return this->m_rhs.evaluate();
    }

    double get_rhs_value(levi::bool_value<false>) {
        return this->m_rhs.evaluate()(0,0);
    }

public:

    DivisionEvaluable(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs)
        : levi::BinaryOperator<product_type, LeftEvaluable, RightEvaluable>(lhs, rhs, lhs.rows(), rhs.cols(), lhs.name() + "/(" + rhs.name() + ")")
    {
        assert(rhs.rows() == 1 && rhs.cols() == 1);
        this->m_info->type = levi::EvaluableType::Division;
        this->m_info->lhs = lhs;
        this->m_info->rhs = rhs;
    }

    virtual ~DivisionEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::row_type>> row(Eigen::Index row) final {
        return this->m_lhs.row(row) / this->m_rhs;
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::col_type>> col(Eigen::Index col) final {
        return this->m_lhs.col(col) / this->m_rhs;
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<product_type>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return this->m_lhs(row, col) / this->m_rhs;
    }

    virtual const product_type& evaluate() final {

        this->m_evaluationBuffer = this->m_lhs.evaluate() / get_rhs_value(levi::bool_value<std::is_arithmetic<typename RightEvaluable::matrix_type>::value>());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<product_type>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {

        levi::ExpressionComponent<typename levi::Evaluable<product_type>::derivative_evaluable> derivative;

        derivative = (this->m_lhs * this->m_rhs.pow(-1)).getColumnDerivative(column, variable);

        return derivative;
    }
};
template <class LeftEvaluable, class RightEvaluable>
levi::DivisionEvaluable<LeftEvaluable, RightEvaluable>::~DivisionEvaluable() { }

/**
 * @brief The SkewEvaluable.
 *
 * It allows to compute the skew symmetric matrix out of a three dimensional vector.
 */
template <typename EvaluableT>
class levi::SkewEvaluable : public levi::UnaryOperator<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>, EvaluableT> {

    Eigen::Matrix<typename EvaluableT::value_type, 3, 1> m_vector;
    levi::ExpressionComponent<levi::TwoElementsMatrix<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>> m_col0, m_col1, m_col2;

    typedef typename EvaluableT::value_type Scalar;
    typedef levi::ExpressionComponent<levi::IdentityEvaluable<Scalar>> ScalarIdentity;

public:

    SkewEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, int)
        : levi::UnaryOperator<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>, EvaluableT>(expression, "skew(" + expression.name() + ")")
          , m_col0(3, 3, levi::Triplet<Scalar>({1, 2, ScalarIdentity()}), levi::Triplet<Scalar>({2, 1, -ScalarIdentity()}), "LeviCivita_ij0")
          , m_col1(3, 3, levi::Triplet<Scalar>({0, 2, -ScalarIdentity()}), levi::Triplet<Scalar>({2, 0, ScalarIdentity()}), "LeviCivita_ij1")
          , m_col2(3, 3, levi::Triplet<Scalar>({0, 1, ScalarIdentity()}), levi::Triplet<Scalar>({1, 0, -ScalarIdentity()}), "LeviCivita_ij2")
    {
        assert((expression.rows() == 3) && (expression.cols() == 1));
    }

    virtual ~SkewEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>::row_type>> row(Eigen::Index row) final {
        return -m_col0.row(row) * this->m_expression(0,0) - m_col1.row(row) * this->m_expression(1,0) - m_col2.row(row) * this->m_expression(2,0);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>::col_type>> col(Eigen::Index col) final {
        return -m_col0.col(col) * this->m_expression(0,0) - m_col1.col(col) * this->m_expression(1,0) - m_col2.col(col) * this->m_expression(2,0);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return -m_col0(row,col) * this->m_expression(0,0) - m_col1(row,col) * this->m_expression(1,0) - m_col2(row,col) * this->m_expression(2,0);
    }

    virtual const Eigen::Matrix<typename EvaluableT::value_type, 3, 3>& evaluate() final {
        m_vector = this->m_expression.evaluate();

        this->m_evaluationBuffer(0,0) = 0.0;
        this->m_evaluationBuffer(0,1) = -m_vector[2];
        this->m_evaluationBuffer(0,2) = m_vector[1];
        this->m_evaluationBuffer(1,0) = m_vector[2];
        this->m_evaluationBuffer(1,1) = 0.0;
        this->m_evaluationBuffer(1,2) = -m_vector[0];
        this->m_evaluationBuffer(2,0) = -m_vector[1];
        this->m_evaluationBuffer(2,1) = m_vector[0];
        this->m_evaluationBuffer(2,2) = 0.0;

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                                                                   std::shared_ptr<levi::VariableBase> variable) final {

        assert( column < 3);

        if (!this->m_expression.isDependentFrom(variable)) {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>::derivative_evaluable::matrix_type>>(this->rows(), this->cols(), "d" + this->name() + "/d" + variable->variableName());
        }

        if (column == 0) {
            return m_col0 * this->m_expression.getColumnDerivative(0, variable);
        }

        if (column == 1) {
            return m_col1 * this->m_expression.getColumnDerivative(0, variable);
        }

        if (column == 2) {
            return m_col2 * this->m_expression.getColumnDerivative(0, variable);
        }

        return levi::ExpressionComponent<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>::derivative_evaluable>();
    }
};
template <typename EvaluableT>
levi::SkewEvaluable<EvaluableT>::~SkewEvaluable() { }

template <typename EvaluableT>
class levi::TransposeEvaluable : public levi::UnaryOperator<typename levi::transpose_type<EvaluableT>::type, EvaluableT> {

    std::vector<levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::value_type>::derivative_evaluable>> m_derivatives;

public:

    TransposeEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, int)
        : levi::UnaryOperator<typename levi::transpose_type<EvaluableT>::type, EvaluableT>(expression, expression.cols(), expression.rows(), expression.name() + "^T")
    {
        this->m_info->type = levi::EvaluableType::Transpose;
        this->m_info->lhs = expression;
    }

    virtual ~TransposeEvaluable() final;

    typedef levi::Evaluable<typename levi::transpose_type<EvaluableT>::type> ThisEvaluable;

    virtual levi::ExpressionComponent<levi::Evaluable<typename ThisEvaluable::row_type>> row(Eigen::Index row) final {
        return this->m_expression.col(row).transpose();
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename ThisEvaluable::col_type>> col(Eigen::Index col) final {
        return this->m_expression.row(col).transpose();
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename ThisEvaluable::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        return this->m_expression(col, row);
    }

    virtual const typename levi::Evaluable<typename levi::transpose_type<EvaluableT>::type>::matrix_type & evaluate() final {
        this->m_evaluationBuffer.lazyAssign(this->m_expression.evaluate().transpose());

        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<typename levi::transpose_type<EvaluableT>::type>::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {
        m_derivatives.resize(this->rows());

        for (Eigen::Index j = 0; j < this->rows(); ++j) {
            m_derivatives[j] = this->m_expression(column, j).getColumnDerivative(0, variable);
        }

        return levi::ComposeByRows(m_derivatives, "d(" + this->name() + ")/d" + variable->variableName());
    }
};
template <typename EvaluableT>
levi::TransposeEvaluable<EvaluableT>::~TransposeEvaluable() { }

template <typename EvaluableT>
class levi::VeeEvaluable : public levi::UnaryOperator<Eigen::Matrix<typename EvaluableT::value_type, 3, 1>, EvaluableT> {

    std::unordered_map<Eigen::Index, std::pair<Eigen::Index, Eigen::Index>> m_indexMap;
    std::vector<levi::ExpressionComponent<typename levi::Evaluable<typename EvaluableT::value_type>::derivative_evaluable>> m_derivatives;


public:

    VeeEvaluable(const levi::ExpressionComponent<EvaluableT>& expression, int)
        : levi::UnaryOperator<Eigen::Matrix<typename EvaluableT::value_type, 3, 1>, EvaluableT>(expression, 3, 1, expression.name() + "^vee")
    {
        this->m_info->lhs = expression;
        m_indexMap[0] = std::make_pair(2,1);
        m_indexMap[1] = std::make_pair(0,2);
        m_indexMap[2] = std::make_pair(1,0);
        m_derivatives.resize(3);
    }

    virtual ~VeeEvaluable() final;

    typedef levi::Evaluable<typename levi::transpose_type<EvaluableT>::type> ThisEvaluable;

    virtual levi::ExpressionComponent<levi::Evaluable<typename ThisEvaluable::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        assert(col == 0);
        levi::unused(col);

        std::pair<Eigen::Index, Eigen::Index> indices = m_indexMap[row];
        return this->m_expression(indices.first, indices.second);
    }

    virtual const Eigen::Matrix<typename EvaluableT::value_type, 3, 1> & evaluate() final {
        const typename EvaluableT::matrix_type original = this->m_expression.evaluate();
        this->m_evaluationBuffer(0) = original(2,1);
        this->m_evaluationBuffer(1) = original(0,2);
        this->m_evaluationBuffer(2) = original(1,0);
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 1>>::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {

        levi::unused(column);

        m_derivatives[0] = this->m_expression(2,1).getColumnDerivative(0, variable);
        m_derivatives[1] = this->m_expression(0,2).getColumnDerivative(0, variable);
        m_derivatives[2] = this->m_expression(1,0).getColumnDerivative(0, variable);

        return levi::ComposeByRows(m_derivatives, "d(" + this->name() + ")/d" + variable->variableName());
    }
};
template <typename EvaluableT>
levi::VeeEvaluable<EvaluableT>::~VeeEvaluable() { }

#endif // LEVI_OPERATORS_H
