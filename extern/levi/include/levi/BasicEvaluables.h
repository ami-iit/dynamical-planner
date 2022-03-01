/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_BASICEVALUABLES_H
#define LEVI_BASICEVALUABLES_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Evaluable.h>
#include <levi/VariableBase.h>
#include <levi/Assignable.h>

template<typename Scalar>
struct levi::TripletStruct{
    Eigen::Index row;
    Eigen::Index col;
    levi::ExpressionComponent<levi::Evaluable<Scalar>> value;
};

/**
 * @brief The ConstantEvaluable
 *
 * Evaluable containing a simple constant matrix
 *
 */
template <typename Matrix>
class levi::ConstantEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> : public levi::Evaluable<Matrix>{
public:

    ConstantEvaluable(std::string name)
        : levi::Evaluable<Matrix>(name)
    {
        this->m_info->type = levi::EvaluableType::Constant;
    }

    ConstantEvaluable(const Matrix& constant, std::string name)
        : levi::Evaluable<Matrix>(constant, name)
    {
        this->m_info->type = levi::EvaluableType::Constant;
    }

    virtual ~ConstantEvaluable() final;

    virtual const Matrix& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Matrix>::derivative_evaluable> getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::derivative_evaluable::matrix_type>>(this->rows(), variable->dimension());
    }

};
template <typename Matrix>
levi::ConstantEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>::~ConstantEvaluable() { }

/**
 * @brief The ConstantEvaluable
 *
 * Evaluable containing a simple constant matrix
 *
 */
template <typename Scalar>
class levi::ConstantEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::Evaluable<Scalar>{
public:

    ConstantEvaluable(const Scalar& constant)
        : levi::Evaluable<Scalar>(constant)
    {
        this->m_info->type = levi::EvaluableType::Constant;
    }

    ConstantEvaluable(const Scalar& constant, const std::string& name)
        : levi::Evaluable<Scalar>(constant, name)
    {
        this->m_info->type = levi::EvaluableType::Constant;
    }

    virtual ~ConstantEvaluable() final;

    virtual const Scalar& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Scalar>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Scalar>::derivative_evaluable::matrix_type>>(1, variable->dimension());
    }

};
template <typename Scalar>
levi::ConstantEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~ConstantEvaluable() { }

/**
 * @brief The MutableEvaluable
 *
 * Evaluable containing a simple matrix, which can be assigned through the operator =
 *
 */
template <typename Matrix>
class levi::MutableEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> : public levi::Assignable<Matrix>
{
public:

    MutableEvaluable(std::string name)
        : levi::Assignable<Matrix>(name)
    { }

    MutableEvaluable(const Matrix& constant, std::string name)
        : levi::Assignable<Matrix>(constant, name)
    { }

    MutableEvaluable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Assignable<Matrix>(rows, cols, name)
    { }

    virtual ~MutableEvaluable() final;

    virtual const Matrix& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Matrix>::derivative_evaluable> getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::derivative_evaluable::matrix_type>>(this->rows(), variable->dimension());
    }
};
template <typename Matrix>
levi::MutableEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>::~MutableEvaluable() { }

/**
 * @brief The ConstantEvaluable
 *
 * Evaluable containing a simple scalar, which can be assigned through the operator =
 *
 */
template <typename Scalar>
class levi::MutableEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::Assignable<Scalar>
{
public:

    MutableEvaluable(const Scalar& constant)
        : levi::Assignable<Scalar>(constant)
    { }

    MutableEvaluable(const Scalar& constant, const std::string& name)
        : levi::Assignable<Scalar>(constant, name)
    { }

    virtual ~MutableEvaluable() final;

    virtual const Scalar& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Scalar>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Scalar>::derivative_evaluable::matrix_type>>(1, variable->dimension());
    }
};
template <typename Scalar>
levi::MutableEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~MutableEvaluable() { }

/**
 * @brief The NullEvaluable
 *
 * Evaluable containing a simple matrix made of zeros.
 *
 */
template <typename Matrix>
class levi::NullEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> : public levi::Evaluable<Matrix>{
public:

    NullEvaluable()
        : levi::Evaluable<Matrix>("0")
    {
        this->m_evaluationBuffer.setZero();
        this->m_info->type = levi::EvaluableType::Null;
    }

    NullEvaluable(Eigen::Index rows, Eigen::Index cols)
        : levi::Evaluable<Matrix>(rows, cols, "0")
    {
        this->m_evaluationBuffer.setZero();
        this->m_info->type = levi::EvaluableType::Null;
    }

    NullEvaluable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Evaluable<Matrix>(rows, cols, name)
    {
        this->m_evaluationBuffer.setZero();
        this->m_info->type = levi::EvaluableType::Null;
    }

    virtual ~NullEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::row_type>> row(Eigen::Index row) final {
        levi::unused(row);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols());
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::col_type>> col(Eigen::Index col) final {
        levi::unused(col);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1);
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::value_type>> element(Eigen::Index row, Eigen::Index col) {
        levi::unused(row, col);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::value_type>>();
    }

    virtual const Matrix& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Matrix>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::derivative_evaluable::matrix_type>>(this->rows(), variable->dimension());
    }
};
template <typename Matrix>
levi::NullEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>::~NullEvaluable() { }

/**
 * @brief The NullEvaluable
 *
 * Evaluable containing a zero.
 *
 */
template <typename Scalar>
class levi::NullEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::Evaluable<Scalar>{
public:

    NullEvaluable()
        : levi::Evaluable<Scalar>(0)
    {
        this->m_info->type = levi::EvaluableType::Null;
    }

    NullEvaluable(Eigen::Index rows, Eigen::Index cols)
        : levi::Evaluable<Scalar>(0)
    {
        this->m_info->type = levi::EvaluableType::Null;
        levi::unused(rows, cols);
        assert(rows == 1 && cols == 1 && "You asked for a null evaluable of type double. It must have only one row and one column.");
    }

    virtual ~NullEvaluable() final;

    virtual const Scalar& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Scalar>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Scalar>::derivative_evaluable::matrix_type>>(1, variable->dimension());
    }

};
template <typename Scalar>
levi::NullEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~NullEvaluable() { }

/**
 * @brief The IdentityEvaluable
 *
 * Evaluable containing a simple identity matrix.
 *
 */
template <typename Matrix>
class levi::IdentityEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> : public levi::Evaluable<Matrix>{
public:

    IdentityEvaluable()
        : levi::Evaluable<Matrix>("I")
    {
        this->m_info->type = levi::EvaluableType::Identity;
        this->m_evaluationBuffer.setIdentity();
    }

    IdentityEvaluable(Eigen::Index rows, Eigen::Index cols)
        : levi::Evaluable<Matrix>(rows, cols, "I")
    {
        this->m_info->type = levi::EvaluableType::Identity;
        this->m_evaluationBuffer.setIdentity();
    }

    IdentityEvaluable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Evaluable<Matrix>(rows, cols, name)
    {
        this->m_info->type = levi::EvaluableType::Identity;
        this->m_evaluationBuffer.setIdentity();
    }

    virtual ~IdentityEvaluable() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::row_type>> row(Eigen::Index row) final {
        levi::Triplet<typename levi::Evaluable<Matrix>::value_type> nonZero = {0, row, levi::ExpressionComponent<levi::IdentityEvaluable<typename levi::Evaluable<Matrix>::value_type>>()};
        return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols(), nonZero,
                                                                                                          "e_" + std::to_string(row) + "^T");
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::col_type>> col(Eigen::Index col) final {
        levi::Triplet<typename levi::Evaluable<Matrix>::value_type> nonZero = {col, 0, levi::ExpressionComponent<levi::IdentityEvaluable<typename levi::Evaluable<Matrix>::value_type>>()};
        return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1, nonZero,
                                                                                                          "e_" + std::to_string(col));

    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::value_type>> element(Eigen::Index row, Eigen::Index col) {
        if (row != col) {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::value_type>>();
        } else {
            return levi::ExpressionComponent<levi::IdentityEvaluable<typename levi::Evaluable<Matrix>::value_type>>();
        }
    }

    virtual const Matrix& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Matrix>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::derivative_evaluable::matrix_type>>(this->rows(), variable->dimension());
    }

};
template <typename Matrix>
levi::IdentityEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>::~IdentityEvaluable() { }

/**
 * @brief The IdentityEvaluable
 *
 * Evaluable containing a 1 of a specified type.
 *
 */
template <typename Scalar>
class levi::IdentityEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::Evaluable<Scalar>{
public:

    IdentityEvaluable()
        : levi::Evaluable<Scalar>(1)
    {
        this->m_info->type = levi::EvaluableType::Identity;
    }

    virtual ~IdentityEvaluable() final;

    virtual const Scalar& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Scalar>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Scalar>::derivative_evaluable::matrix_type>>(1, variable->dimension());
    }

};
template <typename Scalar>
levi::IdentityEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~IdentityEvaluable() { }

template <typename Matrix, typename>
class levi::SingleElementMatrix : public levi::Evaluable<Matrix> {
    levi::Triplet<typename Matrix::value_type> m_element;
public:

    SingleElementMatrix(Eigen::Index rows, Eigen::Index cols, const levi::Triplet<typename Matrix::value_type>& element, const std::string& name)
        : levi::Evaluable<Matrix>(rows, cols, name)
          , m_element(element)
    {
        this->m_info->type = levi::EvaluableType::Constant;

        this->m_evaluationBuffer.setZero();
        auto valueCopy = element.value;
        assert(valueCopy.isValidExpression() && "The specified element is invalid");
        this->m_evaluationBuffer(element.row, element.col) = valueCopy.evaluate();
    }

    virtual ~SingleElementMatrix() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::row_type>> row(Eigen::Index row) final {
        if (row == m_element.row) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> nonZero({0, m_element.col, m_element.value});
            return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols(), nonZero,
                                                                                                              "[" + this->name() +"](" + std::to_string(row) + ",:)");
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols());
        }
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::col_type>> col(Eigen::Index col) final {
        if (col == m_element.col) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> nonZero({m_element.row, 0, m_element.value});
            return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1, nonZero,
                                                                                                              "[" + this->name() +"](:," +
                                                                                                                  std::to_string(col) + ")");
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1);
        }
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        if ((row == m_element.row) && (col == m_element.col)) {
            return m_element.value;
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::value_type>>();
        }
    }

    virtual const Matrix& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Matrix>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::derivative_evaluable::matrix_type>>(this->rows(), variable->dimension());
    }
};
template <typename Matrix, typename empty>
levi::SingleElementMatrix<Matrix, empty>::~SingleElementMatrix() { }

template <typename Matrix, typename>
class levi::TwoElementsMatrix : public levi::Evaluable<Matrix> {
    levi::Triplet<typename Matrix::value_type> m_firstElement;
    levi::Triplet<typename Matrix::value_type> m_secondElement;

public:

    TwoElementsMatrix(Eigen::Index rows, Eigen::Index cols, const levi::Triplet<typename Matrix::value_type>& firstElement,
                      const levi::Triplet<typename Matrix::value_type>& secondElement, const std::string& name)
        : levi::Evaluable<Matrix>(rows, cols, name)
          , m_firstElement(firstElement)
          , m_secondElement(secondElement)
    {
        this->m_info->type = levi::EvaluableType::Constant;

        this->m_evaluationBuffer.setZero();
        auto firstValueCopy = firstElement.value;
        assert(firstValueCopy.isValidExpression() && "The first element is invalid");
        this->m_evaluationBuffer(firstElement.row, firstElement.col) = firstValueCopy.evaluate();
        auto secondValueCopy = secondElement.value;
        assert(secondValueCopy.isValidExpression() && "The second element is invalid");
        this->m_evaluationBuffer(secondElement.row, secondElement.col) = secondValueCopy.evaluate();
    }

    virtual ~TwoElementsMatrix() final;

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::row_type>> row(Eigen::Index row) final {
        if (row == m_firstElement.row && row == m_secondElement.row) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> firstNonZero({0, m_firstElement.col, m_firstElement.value});
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> secondNonZero({0, m_secondElement.col, m_secondElement.value});
            return levi::ExpressionComponent<TwoElementsMatrix<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols(), firstNonZero, secondNonZero,
                                                                                                            "[" + this->name() +"](" + std::to_string(row) + ",:)");
        } else if (row == m_firstElement.row) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> firstNonZero({0, m_firstElement.col, m_firstElement.value});
            return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols(), firstNonZero,
                                                                                                            "[" + this->name() +"](" + std::to_string(row) + ",:)");
        } else if (row == m_secondElement.row) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> secondNonZero({0, m_secondElement.col, m_secondElement.value});
            return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols(), secondNonZero,
                                                                                                              "[" + this->name() +"](" + std::to_string(row) + ",:)");
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::row_type>>(1, this->cols());
        }
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::col_type>> col(Eigen::Index col) final {
        if (col == m_firstElement.col && col == m_secondElement.col) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> firstNonZero({m_firstElement.row, 0, m_firstElement.value});
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> secondNonZero({m_secondElement.row, 0, m_secondElement.value});
            return levi::ExpressionComponent<TwoElementsMatrix<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1, firstNonZero, secondNonZero,
                                                                                                            "[" + this->name() +"](:," +
                                                                                                                std::to_string(col) + ")");
        } else if (col == m_firstElement.col) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> firstNonZero({m_firstElement.row, 0, m_firstElement.value});
            return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1, firstNonZero,
                                                                                                              "[" + this->name() +"](:," +
                                                                                                                  std::to_string(col) + ")");
        } else if (col == m_secondElement.col) {
            levi::Triplet<typename levi::Evaluable<Matrix>::value_type> secondNonZero({m_secondElement.row, 0, m_secondElement.value});
            return levi::ExpressionComponent<SingleElementMatrix<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1, secondNonZero,
                                                                                                              "[" + this->name() +"](:," +
                                                                                                                  std::to_string(col) + ")");
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::col_type>>(this->rows(), 1);
        }
    }

    virtual levi::ExpressionComponent<levi::Evaluable<typename levi::Evaluable<Matrix>::value_type>> element(Eigen::Index row, Eigen::Index col) final {
        if ((row == m_firstElement.row) && (col == m_firstElement.col)) {
            return m_firstElement.value;
        } else if ((row == m_secondElement.row) && (col == m_secondElement.col)) {
            return m_secondElement.value;
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::value_type>>();
        }
    }

    virtual const Matrix& evaluate() final {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Matrix>::derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                                                  std::shared_ptr<levi::VariableBase> variable) final {
        levi::unused(column);
        return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Matrix>::derivative_evaluable::matrix_type>>(this->rows(), variable->dimension());
    }

};
template <typename Matrix, typename empty>
levi::TwoElementsMatrix<Matrix, empty>::~TwoElementsMatrix() { }

#endif // LEVI_BASICEVALUABLES_H
