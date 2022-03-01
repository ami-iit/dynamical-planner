/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_EVALUABLE_H
#define LEVI_EVALUABLE_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/TypeDetector.h>
#include <levi/VariableBase.h>
#include <levi/Registrar.h>

/**
 * @brief Evaluable class (for matrix type)
 *
 * The evaluable class defines a block which can be evaluated. It can be a matrix or a scalar (double, float,..).
 * User can define his own evaluable by publicly inheriting from this class and overload the evaluate() method.
 * It will be on the user to define which value the evaluable will take.
 * The method getColumnDerivative has to be overloaded in order to specify the derivative.
 *
 * Take a look at the ExpressionComponent documentation to place an evaluable into an expression.
 */
template<typename Matrix>
class levi::Evaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> : public levi::Registrar {

    /**
     * @brief Name of the evaluable
     */
    std::string m_name;

    /**
     * @brief Mutex for evaluation
     */
    std::mutex m_evaluationMutex;

public:

    /**
     * @brief The storage type used in the evaluable.
     */
    typedef Matrix matrix_type;

    /**
     * @brief The basic type used into the matrix.
     */
    typedef typename Matrix::value_type value_type;

    /**
     * @brief Number of rows known at compile time (or Eigen::Dynamic, i.e. -1, in case of dynamic types)
     */
    static const Eigen::Index rows_at_compile_time = Matrix::RowsAtCompileTime;

    /**
     * @brief Number of cols known at compile time (or Eigen::Dynamic, i.e. -1, in case of dynamic types)
     */
    static const Eigen::Index cols_at_compile_time = Matrix::ColsAtCompileTime;

    /**
     * @brief The matrix type used to store a single row.
     */
    typedef Eigen::Matrix<value_type, 1, cols_at_compile_time> row_type;

    /**
     * @brief The matrix type used to store a single column.
     */
    typedef Eigen::Matrix<value_type, rows_at_compile_time, 1> col_type;

    /**
     * @brief The Evaluable type used to store the column derivative.
     */
    typedef Evaluable<Eigen::Matrix<value_type, rows_at_compile_time, Eigen::Dynamic>> derivative_evaluable;

    class EvaluableInfo{

        friend class Evaluable;

        EvaluableInfo() { }

    public:
        typedef levi::Evaluable<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> operands_evaluable;
        levi::EvaluableType type;
        levi::BlockType block;
        value_type exponent;
        levi::ExpressionComponent<operands_evaluable> lhs, rhs;
        size_t hash;

        template <typename OtherInfo>
        void copy(const OtherInfo& other) {
            type = other.type;
            block = other.block;
            exponent = other.exponent;
            lhs = other.lhs;
            rhs = other.rhs;
            hash = other.hash;
        }
    };

protected:

    /**
     * @brief The buffer used to store the values of the evaluable.
     */
    Matrix m_evaluationBuffer;

    using DerivativeMap = std::unordered_map<std::string, std::vector<levi::ExpressionComponent<derivative_evaluable>>>;

    using DerivativeMapKey = std::pair<std::string, std::vector<levi::ExpressionComponent<derivative_evaluable>>>;

    /**
     * @brief Map that contains the expressions of derivatives already computed
     */
    DerivativeMap m_derivativeBuffer;

    /**
     * @brief Contains basic infos about the evaluable (useful to expand the evaluation tree)
     */
    EvaluableInfo* m_info;

public:

    Evaluable() = delete;

    /**
     * @brief Constructor
     * @param name Name of the evaluable
     *
     * The evaluation buffer (m_evaluationBuffer) is not initialized nor resized.
     */
    Evaluable(const std::string& name)
        : m_name(name)
        , m_info(new EvaluableInfo)
    {
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    /**
     * @brief Constructor
     * @param rows Number of rows of the evaluable.
     * @param cols Number of columns of the evaluable.
     * @param name Name of the evaluable.
     */
    Evaluable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : m_name(name)
        , m_evaluationBuffer(rows, cols)
        , m_info(new EvaluableInfo)
    {
        m_evaluationBuffer.setZero();
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    /**
     * @brief Constructor
     * @param initialValue Initial value of the evaluation buffer.
     * @param name Name of the evaluable.
     */
    Evaluable(const Matrix& initialValue, const std::string& name)
        : m_name(name)
        , m_evaluationBuffer(initialValue)
        , m_info(new EvaluableInfo)
    {
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    template <typename OtherMatrix>
    Evaluable(const Evaluable<OtherMatrix>& other) = delete;

    template <typename OtherMatrix>
    Evaluable(Evaluable<OtherMatrix>&& other) = delete;

    virtual ~Evaluable();

    /**
     * @brief Number of rows of the evaluable
     * @return the number of rows of the evaluable
     */
    Eigen::Index rows() const {
        return m_evaluationBuffer.rows();
    }

    /**
     * @brief Number of cols of the evaluable
     * @return the number of cols of the evaluable
     */
    Eigen::Index cols() const {
        return m_evaluationBuffer.cols();
    }

    /**
     * @brief Return the name of the evaluable
     * @return The name of the evaluable
     */
    std::string name() const {
        return m_name;
    }

    const EvaluableInfo& info() const {
        return *m_info;
    }

    virtual levi::ExpressionComponent<levi::Evaluable<row_type>> row(Eigen::Index row) {
        levi::unused(row);
        return levi::ExpressionComponent<levi::Evaluable<row_type>>();
    }

    virtual levi::ExpressionComponent<levi::Evaluable<col_type>> col(Eigen::Index col) {
        levi::unused(col);
        return levi::ExpressionComponent<levi::Evaluable<col_type>>();
    }

    virtual levi::ExpressionComponent<levi::Evaluable<value_type>> element(Eigen::Index row, Eigen::Index col) {
        levi::unused(row, col);
        return levi::ExpressionComponent<levi::Evaluable<value_type>>();
    }

    /**
     * @brief Evaluate the evaluable keeping track of the caller
     *
     * @return const reference to the evaluation buffer.
     */
    const Matrix& evaluateID(size_t callerID) {
        std::lock_guard<std::mutex> lock(m_evaluationMutex);

        if (callerID < m_evaluationRegister.size()) {
            if (this->isNew(callerID)) {
                if (this->m_alreadyComputed) {
                    this->m_evaluationRegister[callerID] = this->m_isNewCounter;
                    this->dependenciesChecked();
                    return m_evaluationBuffer;
                } else {
                    const Matrix& output = evaluate();
                    this->m_evaluationRegister[callerID] = this->m_isNewCounter;
                    this->m_alreadyComputed = true;
                    this->dependenciesChecked();
                    return output;
                }
            } else {
                return m_evaluationBuffer;
            }
        }
        return evaluate();
    }

    /**
     * @brief Evaluate the evaluable
     *
     * User should override this method to define custom evaluable.
     *
     * @return const reference to the evaluation buffer.
     */
    virtual const Matrix& evaluate() = 0;

    /**
     * @brief Get the derivative of a specified column with respect to a specified variable
     * @param column The index with respect to the derivative has to be computed.
     * @param variable The variable with respect to the derivative has to be computed.
     * @return An expression containing an evaluable of type derivative_evaluable.
     */
    virtual levi::ExpressionComponent<derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
        levi::unused(column, variable);
        return levi::ExpressionComponent<derivative_evaluable>();
    }

    /**
     * @brief Get the derivative of a specified column with respect to a specified variable
     * @param column The index with respect to the derivative has to be computed.
     * @param variable The variable with respect to the derivative has to be computed.
     * @return An expression containing an evaluable of type derivative_evaluable.
     *
     * @note This will be the method called by ExpressionComponent. This allow to cache the expressions of derivatives which have been already computed.
     * Override this method if you want to avoid caching expressions.
     */
    virtual levi::ExpressionComponent<derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                std::shared_ptr<levi::VariableBase> variable) {
        typename DerivativeMap::iterator element = m_derivativeBuffer.find(variable->variableName());

        levi::ExpressionComponent<derivative_evaluable> emptyExpression;

        if (element == m_derivativeBuffer.end()) {
            DerivativeMapKey newPair;
            newPair.first = variable->variableName();
            newPair.second.resize(this->cols(), emptyExpression);
            auto insertedElement = m_derivativeBuffer.insert(newPair);
            element = insertedElement.first;
        }

        if (!(element->second.at(column).isValidExpression())) {
            if (isDependentFrom(variable)) {
                element->second.at(column) = getNewColumnDerivative(column, variable);
                assert((element->second.at(column).rows() == rows()) && (element->second.at(column).cols() == variable->dimension())
                       && "Wrong dimensions when retrieving new derivative.");
            } else {
                element->second.at(column) = levi::ExpressionComponent<levi::NullEvaluable<typename derivative_evaluable::matrix_type>>(rows(),
                                                                                                                                        variable->dimension(),
                                                                                                                                        "d(" + name() + ")/d" + variable->variableName());
            }
        }

        return element->second.at(column);
    }

    /**
     * @brief Clears the cache of derivatives
     */
    virtual void clearDerivativesCache() {
        m_derivativeBuffer.clear();
    }

};
template<typename Matrix>
levi::Evaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>::~Evaluable()
{
    if (m_info) {
        delete m_info;
        m_info = nullptr;
    }
}

/**
 * @brief Evaluable class (for scalar type)
 *
 * The evaluable class defines a block which can be evaluated. It can be a matrix or a scalar (double, float,..).
 * User can define his own evaluable by publicly inheriting from this class and overload the evaluate() method.
 * It will be on the user to define which value the evaluable will take.
 * The method getColumnDerivative has to be overloaded in order to specify the derivative.
 *
 * Take a look at the ExpressionComponent documentation to place an evaluable into an expression.
 */
template <typename Scalar>
class levi::Evaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::Registrar {

    /**
     * @brief Name of the evaluable
     */
    std::string m_name;

    /**
     * @brief Mutex for evaluation
     */
    std::mutex m_evaluationMutex;

public:

    /**
     * @brief The storage type used in the evaluable.
     *
     * This corresponds also to the value_type, row_type and col_type.
     */
    typedef Scalar matrix_type;

    typedef Scalar value_type;

    static const Eigen::Index rows_at_compile_time = 1;

    static const Eigen::Index cols_at_compile_time = 1;

    typedef Scalar row_type;

    typedef Scalar col_type;

    /**
     * @brief The Evaluable type used to store the column derivative.
     */
    typedef Evaluable<Eigen::Matrix<value_type, 1, Eigen::Dynamic>> derivative_evaluable;

    class EvaluableInfo{

        friend class Evaluable;

        EvaluableInfo() { }

    public:
        typedef levi::Evaluable<Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic>> operands_evaluable;

        levi::EvaluableType type;
        levi::BlockType block;
        value_type exponent;
        levi::ExpressionComponent<operands_evaluable> lhs, rhs;
        size_t hash;

        template <typename OtherInfo>
        void copy(const OtherInfo& other) {
            type = other.type;
            block = other.block;
            exponent = other.exponent;
            lhs = other.lhs;
            rhs = other.rhs;
            hash = other.hash;
        }
    };

protected:

    /**
     * @brief The buffer used to store the values of the evaluable.
     */
    Scalar m_evaluationBuffer;

    using DerivativeMap = std::unordered_map<std::string, std::vector<levi::ExpressionComponent<derivative_evaluable>>>;

    using DerivativeMapKey = std::pair<std::string, std::vector<levi::ExpressionComponent<derivative_evaluable>>>;

    /**
     * @brief Map that contains the expressions of derivatives already computed
     */
    DerivativeMap m_derivativeBuffer;

    EvaluableInfo *m_info;

public:

    Evaluable() = delete;

    /**
     * @brief Constructor
     * @param name Name of the evaluable
     *
     * The evaluation buffer (m_evaluationBuffer) is not initialized.
     */
    Evaluable(const std::string& name)
        : m_name(name)
        , m_info(new EvaluableInfo)
    {
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    /**
     * @brief Constructor (for compatibility with the matrix version)
     * @param rows The number of rows (has to be equal to 1)
     * @param cols The number of cols (has to be equal to 1)
     * @param name Name of the evaluable.
     */
    Evaluable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : m_name(name)
        , m_info(new EvaluableInfo)
    {
        levi::unused(rows, cols);
        assert(rows == 1 && cols == 1);
        m_evaluationBuffer = 0;
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    /**
     * @brief Constructor
     * @param initialValue Initial value of the evaluation buffer.
     * @param name Name of the evaluable.
     */
    Evaluable(const Scalar& initialValue, const std::string& name)
        : m_name(name)
        , m_evaluationBuffer(initialValue)
        , m_info(new EvaluableInfo)
    {
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    /**
     * @brief Constructor
     * @param initialValue Initial value of the evaluation buffer. The name will correspond to this value.
     */
    Evaluable(const Scalar& initialValue)
        : m_name(std::to_string(initialValue))
        , m_evaluationBuffer(initialValue)
        , m_info(new EvaluableInfo)
    {
        m_info->type = levi::EvaluableType::Generic;
        m_info->hash = reinterpret_cast<size_t>(this);
    }

    template <typename OtherMatrix, typename OtherDerivativeEvaluable>
    Evaluable(const Evaluable<OtherMatrix, OtherDerivativeEvaluable>& other) = delete;

    template <typename OtherMatrix, typename OtherDerivativeEvaluable>
    Evaluable(Evaluable<OtherMatrix, OtherDerivativeEvaluable>&& other) = delete;

    virtual ~Evaluable();

    /**
     * @brief Number of rows of the evaluable
     * @return the number of rows of the evaluable
     */
    Eigen::Index rows() const {
        return 1;
    }

    /**
     * @brief Number of cols of the evaluable
     * @return the number of cols of the evaluable
     */
    Eigen::Index cols() const {
        return 1;
    }

    /**
     * @brief Return the name of the evaluable
     * @return The name of the evaluable
     */
    std::string name() const {
        return m_name;
    }

    const EvaluableInfo& info() const {
        return *m_info;
    }

    virtual levi::ExpressionComponent<levi::Evaluable<row_type>> row(Eigen::Index row) {
        levi::unused(row);
        return levi::ExpressionComponent<levi::Evaluable<row_type>>();
    }

    virtual levi::ExpressionComponent<levi::Evaluable<col_type>> col(Eigen::Index col) {
        levi::unused(col);
        return levi::ExpressionComponent<levi::Evaluable<col_type>>();
    }

    virtual levi::ExpressionComponent<levi::Evaluable<value_type>> element(Eigen::Index row, Eigen::Index col) {
        levi::unused(row, col);
        return levi::ExpressionComponent<levi::Evaluable<value_type>>();
    }

    /**
     * @brief Evaluate the evaluable keeping track of the caller
     *
     * @return const reference to the evaluation buffer.
     */
    const Scalar& evaluateID(size_t callerID) {
        std::lock_guard<std::mutex> lock(m_evaluationMutex);

        if (callerID < m_evaluationRegister.size()) {
            if (this->isNew(callerID)) {
                if (this->m_alreadyComputed) {
                    this->m_evaluationRegister[callerID] = this->m_isNewCounter;
                    this->dependenciesChecked();
                    return m_evaluationBuffer;
                } else {
                    const Scalar& output = evaluate();
                    this->m_evaluationRegister[callerID] = this->m_isNewCounter;
                    this->m_alreadyComputed = true;
                    this->dependenciesChecked();
                    return output;
                }
            } else {
                return m_evaluationBuffer;
            }
        }

        return evaluate();
    }

    /**
     * @brief Evaluate the evaluable
     *
     * User should override this method to define custom evaluable.
     *
     * @return const reference to the evaluation buffer.
     */
    virtual const Scalar& evaluate() = 0;

    /**
     * @brief Get the derivative of a specified column with respect to a specified variable
     * @param column The index with respect to the derivative has to be computed.
     * @param variable The variable with respect to the derivative has to be computed.
     * @return An expression containing an evaluable of type derivative_evaluable.
     */
    virtual levi::ExpressionComponent<derivative_evaluable> getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) {
        levi::unused(column, variable);
        return levi::ExpressionComponent<derivative_evaluable>();
    }

    /**
     * @brief Get the derivative of a specified column with respect to a specified variable
     * @param column The index with respect to the derivative has to be computed.
     * @param variable The variable with respect to the derivative has to be computed.
     * @return An expression containing an evaluable of type derivative_evaluable.
     *
     * @note This will be the method called by ExpressionComponent. This allow to cache the expressions of derivatives which have been already computed.
     * Override this method if you want to avoid caching expressions.
     */
    virtual levi::ExpressionComponent<derivative_evaluable> getColumnDerivative(Eigen::Index column,
                                                                                std::shared_ptr<levi::VariableBase> variable) {
        typename DerivativeMap::iterator element = m_derivativeBuffer.find(variable->variableName());

        levi::ExpressionComponent<derivative_evaluable> emptyExpression;

        if (element == m_derivativeBuffer.end()) {
            DerivativeMapKey newPair;
            newPair.first = variable->variableName();
            newPair.second.resize(this->cols(), emptyExpression);
            auto insertedElement = m_derivativeBuffer.insert(newPair);
            element = insertedElement.first;
        }

        if (!(element->second.at(column).isValidExpression())) {
            if (isDependentFrom(variable)) {
                element->second.at(column) = getNewColumnDerivative(column, variable);
                assert((element->second.at(column).rows() == rows()) && (element->second.at(column).cols() == variable->dimension())
                       && "Wrong dimensions when retrieving new derivative.");
            } else {
                element->second.at(column) = levi::ExpressionComponent<levi::NullEvaluable<typename derivative_evaluable::matrix_type>>(rows(),
                                                                                                                                        variable->dimension(),
                                                                                                                                        "d(" + name() + ")/d" + variable->variableName());
            }
        }

        return element->second.at(column);
    }

    /**
     * @brief Clears the cache of derivatives
     */
    virtual void clearDerivativesCache() {
        m_derivativeBuffer.clear();
    }

};
template <typename Scalar>
levi::Evaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~Evaluable()
{
    if (m_info) {
        delete m_info;
        m_info = nullptr;
    }
}

#endif // LEVI_EVALUABLE_H
