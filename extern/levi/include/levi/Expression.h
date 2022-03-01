/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_EXPRESSION_COMPONENT_H
#define LEVI_EXPRESSION_COMPONENT_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>

namespace levi {

    template<bool value, typename T>
    static ExpressionComponent<ConstantEvaluable<T>> build_constant(bool_value<value>, const T& rhs);

    /**
    * @brief Return an expression containing a constant evaluable (arithmetic type)
    */
    template<typename T>
    static ExpressionComponent<ConstantEvaluable<T>> build_constant(bool_value<true>, const T& rhs);

    /**
    * @brief Return an expression containing a constant evaluable (matrix type)
    *
    * The name of the constant will be "UnnamedMatrix"
    */
    template<typename T>
    static ExpressionComponent<ConstantEvaluable<T>> build_constant(bool_value<false>, const T& rhs);

    /**
    * @brief Create an epression from a vector of rows
    * @param rows The vector containing the expressionion composing the rows of the new expression
    * @param name The name of the new expression
    * @return An expression made of the specified rows
    */
    template<int rowsNumber = -1, typename EvaluableT>
    ExpressionComponent<typename levi::ConstructorByRows<EvaluableT, rowsNumber>::composite_evaluable>
    ComposeByRows(const std::vector<ExpressionComponent<EvaluableT>>& rows, const std::string& name);

    /**
     * @brief Create an epression from a vector of cols
     * @param cols The vector containing the expressionion composing the columns of the new expression
     * @param name The name of the new expression
     * @return An expression made of the specified columns
     */
    template<int colsNumber = -1, typename EvaluableT>
    ExpressionComponent<typename levi::ConstructorByCols<EvaluableT, colsNumber>::composite_evaluable>
    ComposeByCols(const std::vector<levi::ExpressionComponent<EvaluableT>>& cols, const std::string& name);

    template <typename LhsT, typename RhsT>
    using ProductOutputEvaluable = levi::Evaluable<typename levi::matrix_product_return<typename LhsT::matrix_type, typename RhsT::matrix_type>::type>;

    template <typename EvaluableLhs, typename EvaluableRhs, typename Enabler = void>
    struct are_commutable : std::false_type {};

    template <typename EvaluableLhs, typename EvaluableRhs>
    struct are_commutable<EvaluableLhs, EvaluableRhs, typename std::enable_if<levi::is_valid_product<EvaluableRhs::rows_at_compile_time,
                                                                                                     EvaluableRhs::cols_at_compile_time,
                                                                                                     EvaluableLhs::rows_at_compile_time,
                                                                                                     EvaluableLhs::cols_at_compile_time>::value>::type> : std::true_type {};

    template <typename Matrix>
    using CompiledElement = std::pair<std::string, levi::ExpressionComponent<levi::Evaluable<Matrix>>>;

    typedef CompiledElement<LEVI_DEFAULT_MATRIX_TYPE> DefaultCompiledElement;

    template <typename Matrix>
    using MultipleExpressionsMap = std::unordered_map<std::string, levi::ExpressionComponent<levi::Evaluable<Matrix>>>;

    typedef MultipleExpressionsMap<LEVI_DEFAULT_MATRIX_TYPE> DefaultMultipleExpressionsMap;

    template <typename Matrix>
    using MultipleCompiledOutputPointer = std::unique_ptr<levi::MultipleCompiledExpressions<levi::Evaluable<Matrix>>>;

    typedef MultipleCompiledOutputPointer<LEVI_DEFAULT_MATRIX_TYPE> DefaultMultipleCompiledOutputPointer;

    template <typename Matrix>
    using MultipleSqueezedOutputPointer = std::unique_ptr<levi::MultipleSqueezedExpressions<levi::Evaluable<Matrix>>>;

    typedef MultipleSqueezedOutputPointer<LEVI_DEFAULT_MATRIX_TYPE> DefaultMultipleSqueezedOutputPointer;

    template <typename Matrix>
    using MultipleExpressionsOutputMap = std::unordered_map<std::string, Matrix>;

    typedef std::unordered_map<std::string, LEVI_DEFAULT_MATRIX_TYPE> DefaultMultipleExpressionsOutputMap;

    template <typename Matrix>
    MultipleCompiledOutputPointer<Matrix> CompileMultipleExpressions(const MultipleExpressionsMap<Matrix>& elements, const std::string& name);

    template <typename Matrix>
    MultipleSqueezedOutputPointer<Matrix> SqueezeMultipleExpressions(const MultipleExpressionsMap<Matrix>& elements);

    template<typename LeftEvaluable, typename RightEvaluable>
    ExpressionComponent<Evaluable<Eigen::Matrix<typename levi::scalar_product_return<typename LeftEvaluable::value_type,
                                                                                     typename RightEvaluable::value_type>::type,
                                                LeftEvaluable::rows_at_compile_time, Eigen::Dynamic>>>
    MatrixProductDerivativeExpression(const levi::ExpressionComponent<LeftEvaluable>& lhs, const levi::ExpressionComponent<RightEvaluable>& rhs, std::shared_ptr<levi::VariableBase> variable);
}

/**
 *@brief The ExpressionComponent class
 *
 * An ExpressionComponent stores a pointer to an evaluable.
 * This allows to perform operations between evaluable, without carrying around the evaluation buffers.
 */
template <class EvaluableT>
class levi::ExpressionComponent {

    template <class OtherEvaluable>
    friend class ExpressionComponent;

    /**
     * @brief Shared pointer to the evaluable.
     */
    std::shared_ptr<EvaluableT> m_evaluable;

    /**
     * @brief CallerID associated to this expression for the specific evaluable;
     */
    size_t m_callerID;

    /**
     * Template declaration for a helper method used when calling the default constructor of the ExpressionComponent
     */
    template<bool value>
    void default_constructor(levi::bool_value<value>);

    /**
     * Template specialization for the case in which the default constructor is available.
     */
    void default_constructor(levi::bool_value<true>);

    /**
     * Template specialization for the case in which the default constructor is *not* available.
     */
    void default_constructor(levi::bool_value<false>);

    /**
     * Template declaration of a function used to perform casting when during assignements.
     */
    template<bool value, typename OtherEvaluable>
    void casted_assignement(levi::bool_value<value>, const ExpressionComponent<OtherEvaluable>& other);

    /**
     * Template specialization in case the two pointers are directly castable.
     */
    template<typename OtherEvaluable>
    void casted_assignement(levi::bool_value<true>, const levi::ExpressionComponent<OtherEvaluable>& other);

    /**
     * Template specialization in case the two pointers are not directly castable.
     */
    template<typename OtherEvaluable>
    void casted_assignement(levi::bool_value<false>, const levi::ExpressionComponent<OtherEvaluable>& other);

    /**
     * Helper class to check if it is suitable to return *this
     */
    template<bool isSuitable, typename EvaluableOut>
    levi::ExpressionComponent<EvaluableOut> return_this(levi::bool_value<isSuitable>) const;

    template<typename EvaluableOut>
    levi::ExpressionComponent<EvaluableOut> return_this(levi::bool_value<true>) const;

    template<typename EvaluableOut>
    levi::ExpressionComponent<EvaluableOut> return_this(levi::bool_value<false>) const;

    /**
     * Helper class to check if it is suitable to return the rhs
     */
    template<bool isSuitable, typename EvaluableOut, typename EvaluableRhs>
    levi::ExpressionComponent<EvaluableOut> return_rhs(levi::bool_value<isSuitable>, const levi::ExpressionComponent<EvaluableRhs>& rhs) const;

    template<typename EvaluableOut, typename EvaluableRhs>
    levi::ExpressionComponent<EvaluableOut> return_rhs(levi::bool_value<true>, const levi::ExpressionComponent<EvaluableRhs>& rhs) const;

    template<typename EvaluableOut, typename EvaluableRhs>
    levi::ExpressionComponent<EvaluableOut> return_rhs(levi::bool_value<false>, const levi::ExpressionComponent<EvaluableRhs>& rhs) const;

    template<bool value, typename EvaluableRhs>
    levi::ExpressionComponent<levi::ProductOutputEvaluable<EvaluableT, EvaluableRhs>> commute_factors(levi::bool_value<value>, const levi::ExpressionComponent<EvaluableRhs>& rhs) const;

    template<typename EvaluableRhs>
    levi::ExpressionComponent<levi::ProductOutputEvaluable<EvaluableT, EvaluableRhs>> commute_factors(levi::bool_value<true>, const levi::ExpressionComponent<EvaluableRhs>& rhs) const;

    template<typename EvaluableRhs>
    levi::ExpressionComponent<levi::ProductOutputEvaluable<EvaluableT, EvaluableRhs>> commute_factors(levi::bool_value<false>, const levi::ExpressionComponent<EvaluableRhs>& rhs) const;


    /**
     * @brief Check if the current evaluable is null
     */
    bool m_isNull;

    /**
     * @brief Check if the current evaluable is a squared identity
     */
    bool m_isIdentity;

public:

    /**
     * @brief Default Constructor
     *
     * If the evaluable EvaluableT is not default constructible, then the shared pointer is equal to nullptr.
     */
    ExpressionComponent();

    /**
     * @brief Copy constructor from other expressions
     */
    template<class EvaluableOther, typename = typename std::enable_if<!std::is_same<EvaluableT, EvaluableOther>::value>::type>
    ExpressionComponent(const ExpressionComponent<EvaluableOther>& other);

    /**
     * @brief Copy constructor
     */
    ExpressionComponent(const ExpressionComponent<EvaluableT>& other);

    /**
     * @brief Move constructor from other expressions
     */
    template<class EvaluableOther, typename = typename std::enable_if<!std::is_same<EvaluableT, EvaluableOther>::value>::type>
    ExpressionComponent(ExpressionComponent<EvaluableOther>&& other);

    /**
     * @brief Move constructor
     */
    ExpressionComponent(ExpressionComponent<EvaluableT>&& other);

    /**
     * @brief Constructor
     *
     * Use this constructor to instantiate all the custom evaluables.
     * All the input arguments will be passed to the constructor of the Evaluable pointed inside the ExpressionComponent.
     */
    template<class... Args, typename = typename std::enable_if<std::is_constructible<EvaluableT, Args...>::value>::type>
    ExpressionComponent(Args&&... args);

    ~ExpressionComponent();

    /**
     * @brief Weak pointer to the evaluable
     *
     * Useful in case the evaluable has to be accessed after the allocation.
     *
     * @return a weak pointer to the evaluable
     */
    std::weak_ptr<EvaluableT> evaluable() const ;

    const typename EvaluableT::EvaluableInfo &info() const;

    /**
     * @brief Name of the expression
     *
     * Corresponds to the name of the pointed evaluable.
     *
     * @return The name of the expression
     */
    std::string name() const;

    /**
     * @brief Rows of the expression
     *
     * Corresponds to the number of rows of the pointed evaluable.
     *
     * @return The number of rows of the expression
     */
    Eigen::Index rows() const;

    /**
     * @brief Cols of the expression
     *
     * Corresponds to the number of cols of the pointed evaluable.
     *
     * @return The number of cols of the expression
     */
    Eigen::Index cols() const;

    /**
     * @brief Check if the expression has a new value
     * @return True if new
     */
    bool isNew() const;

    /**
     * @brief Check if the expression is null
     * @return True if null
     */
    bool isNull() const;

    /**
     * @brief Evaluate the pointed evaluable.
     *
     * @warning An assert is used to check whether the evaluable can be evaluated. Test the code in debug mode first to avoid segfaults.
     *
     * @return A const reference to the evaluation buffer of the pointed evaluable.
     */
    const typename EvaluableT::matrix_type& evaluate();

    /**
     * @brief Operator +
     * @return An expression which points to an evaluable performing the additions.
     */
    template<class EvaluableRhs>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_sum_return<typename EvaluableT::matrix_type, typename EvaluableRhs::matrix_type>::type>> operator+(const ExpressionComponent<EvaluableRhs>& rhs) const;

    /**
     * @brief Operator +
     *
     * The other addend will be inserted in a ConstantEvaluable containing the rhs.
     *
     * @return An expression which points to an evaluable performing the addition.
     */
    template <typename Matrix>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_sum_return<typename EvaluableT::matrix_type, Matrix>::type>> operator+(const Matrix& rhs) const;

    /**
     * @brief Operator +
     */
    ExpressionComponent<levi::Evaluable<typename EvaluableT::matrix_type>> operator+() const;

    /**
     * @brief Operator -
     *
     * @return An expression which points to an evaluable performing the subtraction.
     */
    template<class EvaluableRhs>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_sum_return<typename EvaluableT::matrix_type, typename EvaluableRhs::matrix_type>::type>> operator-(const ExpressionComponent<EvaluableRhs>& rhs) const;

    /**
     * @brief Operator -
     *
     * The other addend will be inserted in a ConstantEvaluable containing the rhs.
     *
     * @return An expression which points to an evaluable performing the subtraction.
     */
    template <typename Matrix>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_sum_return<typename EvaluableT::matrix_type, Matrix>::type>> operator-(const Matrix& rhs) const;

    /**
     * @brief Operator -
     *
     * @return An expression which points to an evaluable inverting the sign.
     */
    ExpressionComponent<levi::Evaluable<typename EvaluableT::matrix_type>> operator-() const;

    /**
     * @brief Operator *
     *
     * @return An expression which points to an evaluable performing the multiplication.
     */
    template<class EvaluableRhs>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_product_return<typename EvaluableT::matrix_type, typename EvaluableRhs::matrix_type>::type>> operator*(const ExpressionComponent<EvaluableRhs>& rhs) const;

    /**
     * @brief Operator *
     *
     * The other addend will be inserted in a ConstantEvaluable containing the rhs.
     *
     * @return An expression which points to an evaluable performing the multiplication.
     */
    template <typename Matrix>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_product_return<typename EvaluableT::matrix_type, Matrix>::type>> operator*(const Matrix& rhs) const;

    /**
     * @brief Operator /
     *
     * @return An expression which points to an evaluable performing the division by a scalar.
     */
    template<class EvaluableRhs>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_product_return<typename EvaluableT::matrix_type, typename EvaluableRhs::value_type>::type>> operator/(const ExpressionComponent<EvaluableRhs>& rhs) const;

    /**
     * @brief Operator /
     *
     * The other addend will be inserted in a ConstantEvaluable containing the rhs.
     *
     * @return An expression which points to an evaluable performing the division.
     */
    template <typename Scalar>
    ExpressionComponent<levi::Evaluable<typename levi::matrix_product_return<typename EvaluableT::matrix_type, Scalar>::type>> operator/(const Scalar& rhs) const;

    /**
     * @brief Operator ==
     *
     * Checks if the pointed evaluables are the same.
     */
    template <typename EvaluableRhs>
    bool operator==(const levi::ExpressionComponent<EvaluableRhs>& other) const;

    /**
     * @brief Operator !=
     *
     * Checks if the pointed evaluables are different.
     */
    template <typename EvaluableRhs>
    bool operator!=(const levi::ExpressionComponent<EvaluableRhs>& other) const;

    /**
     * @brief Computes the power of the current evaluable
     * @param exponent The exponent of the power
     * @return An expression computing the power
     *
     * @note This can be used only with scalars or 1x1 matrices.
     */
    ExpressionComponent<levi::Evaluable<typename EvaluableT::value_type>> pow(typename EvaluableT::value_type exponent) const;

    /**
     * @brief Assignement operator from other expressions
     *
     * Assigns the current expression to the rhs. If the EvaluableT is not a base class for EvaluableRhs,
     * a new evaluable will be created which will cast the evaluation buffers.
     *
     * @note This changes the pointer and does not affect the pointed evaluable. If no other expression points to the previous evaluable,
     * it will be deleted.
     */
    template<class EvaluableRhs, typename = typename std::enable_if<!std::is_same<EvaluableT, EvaluableRhs>::value>::type>
    void operator=(const ExpressionComponent<EvaluableRhs>& rhs);

    /**
     * @brief Assignement operator
     *
     * @note This changes the pointer and does not affect the pointed evaluable. If no other expression points to the previous evaluable,
     * it will be deleted.
     */
    void operator=(const ExpressionComponent<EvaluableT>& rhs);

    /**
     * @brief Move assignement operator from other expressions
     *
     * Assigns the current expression to the rhs. If the EvaluableT is not a base class for EvaluableRhs,
     * a new evaluable will be created which will cast the evaluation buffers.
     */
    template<class EvaluableRhs, typename = typename std::enable_if<!std::is_same<EvaluableT, EvaluableRhs>::value>::type>
    void operator=(const ExpressionComponent<EvaluableRhs>&& rhs);

    /**
     * @brief Move assignement operator
     */
    void operator=(const ExpressionComponent<EvaluableT>&& rhs);

    /**
     * @brief Assignement operator to a matrix
     *
     * Assigns the evaluable pointed by the current expression to the rhs. This method can be used only if the operator= is defined
     * in the pointed evaluable for the Matrix type.
     */
    template<typename Matrix>
    void operator=(const Matrix& rhs);

    /**
     * @brief Returns an expression corresponding to a row of the current evaluable
     * @param row The index of the row of interest
     * @return An expression pointing to a row of the current evaluable.
     *
     * @Note This has to be considered as read only accessor.
     */
    ExpressionComponent<levi::Evaluable<typename EvaluableT::row_type>> row(Eigen::Index row) const;

    /**
     * @brief Returns an expression corresponding to a column of the current evaluable
     * @param col The index of the column of interest
     * @return An expression pointing to a column of the current evaluable.
     *
     * @Note This has to be considered as read only accessor.
     */
    ExpressionComponent<levi::Evaluable<typename EvaluableT::col_type>> col(Eigen::Index col) const;

    /**
     * @brief Accessor to element (read only)
     * @param row The index of the row of interest
     * @param col The index of the column of interest
     * @return An expression pointing to an element of the current evaluable.
     *
     * @Note This has to be considered as read only accessor.
     *
     */
    ExpressionComponent<levi::Evaluable<typename EvaluableT::value_type>> operator()(Eigen::Index row, Eigen::Index col) const;

    /**
     * @brief Accessor to element (read only)
     * @param startRow The starting row to be considered in the block
     * @param startCol The starting column to be considered in the block
     * @param numberOfRows The number of rows of the block
     * @param numberOfCols The number of columns of the block
     * @return An expression pointing to an evaluable taking only the specified block out of the current evaluable.
     */
    ExpressionComponent<levi::Evaluable<typename levi::dynamic_block_return<typename EvaluableT::matrix_type>::type>> block(Eigen::Index startRow, Eigen::Index startCol, Eigen::Index numberOfRows, Eigen::Index numberOfCols) const;

    /**
     * @brief Accessor to element (read only)
     * @param startRow The starting row to be considered in the block
     * @param startCol The starting column to be considered in the block
     * @return An expression pointing to an evaluable taking only the specified block out of the current evaluable.
     */
    template<unsigned int numberOfRows, unsigned int numberOfCols>
    ExpressionComponent<levi::Evaluable<typename levi::fixed_block_return<typename EvaluableT::matrix_type, numberOfRows, numberOfCols>::type>> block(Eigen::Index startRow, Eigen::Index startCol) const;


    /**
     * @brief Computes the skew symmetric matrix out of the current evaluable
     * @return An expression pointing to an evaluable which computes the skew symmetric matrix
     * @note This works only with three dimensional vectors.
     */
    ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 3>>> skew() const;

    /**
     * @brief Extract the vector defining a skew symmetric matrix
     * @return An expression pointing to an evaluable which performs the vee operator
     * @note This works only with three dimensional squared matrices.
     */
    ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, 3, 1>>> vee() const;


    /**
     * @brief Computes the transpose of the current expression
     * @return An expression whose evaluable computes the transpose of the current evaluable.
     *
     * @note The column derivative needs to call the derivative of every column. This may be expensive.
     */
    ExpressionComponent<levi::Evaluable<typename levi::transpose_type<EvaluableT>::type>> transpose() const;

    /**
     * @brief get a variable from the current expression
     * @return A variable obtained by interpreting the current expression as a variable.
     *
     * @note The expression must have a single column.
     */
    ExpressionComponent<levi::EvaluableVariable<typename EvaluableT::col_type>> asVariable() const;

    /**
     * @brief get a variable from the current expression
     * @return A variable obtained by interpreting the current expression as a variable. This is like a completely
     * new variable, thus the derivative wrt another variable will be zero. Nevertheless, the evaluation will still
     * evaluate the underneath expression.
     *
     * @note The expression must have a single column.
     */
    ExpressionComponent<levi::EvaluableVariable<typename EvaluableT::col_type>> asIndependentVariable() const;

    /**
     * @brief Generates a new expression condensing all the nodes in one. This expression canno be derived.
     * @return A new expression containing the condensed version of the current expression
     */
    ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic>>> squeeze(const std::string &name) const;

    /**
     * @brief Generates a new expression condensing all the nodes in one. This expression canno be derived.
     * @return A new expression containing the condensed version of the current expression
     */
    ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic>>> compile(const std::string &name) const;

    /**
     * @brief Retrieve the column derivative with respect to the specified variable
     *
     * @param column The column of interest
     * @param variable Expression pointing to the variable of interest
     * @return An expression pointing to the column derivative.
     */
    template<typename VariableType>
    ExpressionComponent<typename EvaluableT::derivative_evaluable> getColumnDerivative(Eigen::Index column, const ExpressionComponent<levi::EvaluableVariable<VariableType>>& variable) const;

    /**
     * @brief Retrieve the column derivative with respect to the specified variable
     *
     * @param column The column of interest
     * @param variable Shared pointer to the variable of interest
     * @return An expression pointing to the column derivative.
     */
    ExpressionComponent<typename EvaluableT::derivative_evaluable> getColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) const;

    /**
     * @brief Clears the cache of derivatives
     */
    void clearDerivativesCache();

    /**
     * @brief Check whether the pointed evaluable depends on a specified variable
     * @param variable The variable of interest
     * @return True if dependent
     */
    template<typename VariableType>
    bool isDependentFrom(const ExpressionComponent<levi::EvaluableVariable<VariableType>>& variable) const;

    /**
     * @brief Check whether the pointed evaluable depends on a specified variable
     * @param variable The variable of interest
     * @return True if dependent
     */
    bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) const;

    std::vector<Registrar*> getDependencies() const;

    /**
     * @brief Check if the expression is valid
     * @return True if the expression is valid
     */
    bool isValidExpression() const;

    /**
     * @brief Create an epression from a vector of rows
     * @param rows The vector containing the expressionion composing the rows of the new expression
     * @param name The name of the new expression
     * @return An expression made of the specified rows
     */
    static ExpressionComponent<EvaluableT> ComposeByRows(const std::vector<levi::ExpressionComponent<levi::Evaluable<typename EvaluableT::row_type>>>& rows, std::string name);

    /**
     * @brief Create an epression from a vector of cols
     * @param cols The vector containing the expressionion composing the columns of the new expression
     * @param name The name of the new expression
     * @return An expression made of the specified columns
     */
    static ExpressionComponent<EvaluableT> ComposeByCols(const std::vector<levi::ExpressionComponent<levi::Evaluable<typename EvaluableT::col_type>>>& cols, std::string name);

    template <typename LeftEvaluable, typename RightEvaluable>
    static ExpressionComponent<EvaluableT> Horzcat(const ExpressionComponent<LeftEvaluable>& lhs,
                                                   const ExpressionComponent<RightEvaluable>& rhs,
                                                   const std::string& name);

    template <typename TopEvaluable, typename BottomEvaluable>
    static ExpressionComponent<EvaluableT> Vertcat(const ExpressionComponent<TopEvaluable>& top,
                                                   const ExpressionComponent<BottomEvaluable>& bottom,
                                                   const std::string& name);
};


/**
 * @brief Operator +
 *
 * The other addend will be inserted in a ConstantEvaluable containing the lhs.
 *
 * @return An expression which points to an evaluable performing the addition.
 */
template <typename Matrix, class EvaluableT>
levi::ExpressionComponent<levi::Evaluable<typename levi::matrix_sum_return<typename EvaluableT::matrix_type, Matrix>::type>> operator+(const Matrix& lhs, const levi::ExpressionComponent<EvaluableT> &rhs);

/**
 * @brief Operator -
 *
 * The other addend will be inserted in a ConstantEvaluable containing the lhs.
 *
 * @return An expression which points to an evaluable performing the subrraction.
 */
template <typename Matrix, class EvaluableT>
levi::ExpressionComponent<levi::Evaluable<typename levi::matrix_sum_return<typename EvaluableT::matrix_type, Matrix>::type>> operator-(const Matrix& lhs, const levi::ExpressionComponent<EvaluableT> &rhs);

/**
 * @brief Operator *
 *
 * The other addend will be inserted in a ConstantEvaluable containing the lhs.
 *
 * @return An expression which points to an evaluable performing the multiplication.
 */
template <typename Matrix, class EvaluableT>
levi::ExpressionComponent<levi::Evaluable<typename levi::matrix_product_return<Matrix, typename EvaluableT::matrix_type>::type>> operator*(const Matrix& lhs, const levi::ExpressionComponent<EvaluableT> &rhs);

/**
 * @brief Operator /
 *
 * The other addend will be inserted in a ConstantEvaluable containing the lhs.
 *
 * @return An expression which points to an evaluable performing the division.
 */
template <typename Matrix, class EvaluableT>
levi::ExpressionComponent<levi::Evaluable<typename levi::matrix_product_return<Matrix, typename EvaluableT::value_type>::type>> operator/(const Matrix& lhs, const levi::ExpressionComponent<EvaluableT> &rhs);


#endif // LEVI_EXPRESSION_COMPONENT_H
