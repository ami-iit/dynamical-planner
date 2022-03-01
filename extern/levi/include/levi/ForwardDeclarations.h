/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_FORWARDDECLARATIONS_H
#define LEVI_FORWARDDECLARATIONS_H

#ifndef LEVI_DEFAULT_MATRIX_TYPE
#define LEVI_DEFAULT_MATRIX_TYPE levi::INVALID_TYPE
#endif

#ifndef LEVI_DEFAULT_MATRIX_FIX_TYPE
#define LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, cols) INVALID_FIXED_TYPE<rows, cols>
#endif

#ifndef LEVI_DEFAULT_VECTOR_TYPE
#define LEVI_DEFAULT_VECTOR_TYPE levi::INVALID_TYPE
#endif

#include <type_traits>
#include <cstddef>


namespace levi {

    class INVALID_TYPE;

    template<int rows, int cols>
    class INVALID_FIXED_TYPE;

    class Registrar;

    template<typename Matrix, class Enabler = void>
    class Evaluable;

    template <typename Matrix>
    class Evaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template <typename Scalar>
    class Evaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template <typename Matrix, class Enabler = void>
    class Assignable;

    template <typename Matrix>
    class Assignable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template <typename Scalar>
    class Assignable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    class VariableBase;

    template <typename Vector, class Enabler = void>
    class EvaluableVariable;

    template <typename Vector>
    class EvaluableVariable<Vector, typename std::enable_if<!std::is_arithmetic<Vector>::value>::type>;

    template <typename Scalar>
    class EvaluableVariable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template <class Evaluable>
    class ExpressionComponent;

    template <typename Matrix, class Enabler = void>
    class ConstantEvaluable;

    template <typename Matrix>
    class ConstantEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template <typename Scalar>
    class ConstantEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template <typename Matrix, class Enabler = void>
    class MutableEvaluable;

    template <typename Matrix>
    class MutableEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template <typename Scalar>
    class MutableEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template <typename Matrix, class Enabler = void>
    class NullEvaluable;

    template <typename Matrix>
    class NullEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template <typename Scalar>
    class NullEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template <typename Matrix, class Enabler = void>
    class IdentityEvaluable;

    template <typename Matrix>
    class IdentityEvaluable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template <typename Scalar>
    class IdentityEvaluable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template <typename Matrix, typename = typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>
    class SingleElementMatrix;

    template <typename Matrix, typename = typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>
    class TwoElementsMatrix;

    template <typename MatrixType, typename EvaluableT>
    class UnaryOperator;

    template <typename MatrixType, class LeftEvaluable, class RightEvaluable>
    class BinaryOperator;

    template <class LeftEvaluable, class RightEvaluable>
    class SumEvaluable;

    template <class LeftEvaluable, class RightEvaluable>
    class SubtractionEvaluable;

    template <class EvaluableT>
    class SignInvertedEvaluable;

    template <class LeftEvaluable, class RightEvaluable>
    class ProductEvaluable;

    template <class LeftEvaluable, class RightEvaluable>
    class MatrixProductDerivative;

    template <class EvaluableT>
    class PowEvaluable;

    template <class LeftEvaluable, class RightEvaluable>
    class DivisionEvaluable;

    template<class Evaluable, class Enabler = void>
    class RowEvaluable;

    template <typename EvaluableT>
    class RowEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template <typename EvaluableT>
    class RowEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template<class Evaluable, class Enabler = void>
    class ColEvaluable;

    template <typename EvaluableT>
    class ColEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template <typename EvaluableT>
    class ColEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template<class EvaluableT, class Enabler = void>
    class ElementEvaluable;

    template <typename EvaluableT>
    class ElementEvaluable<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template <typename EvaluableT>
    class ElementEvaluable<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template<class InputEvaluable, class OutputEvaluable, class Enabler = void>
    class BlockEvaluable;

    template <class InputEvaluable, class OutputEvaluable>
    class BlockEvaluable<InputEvaluable, OutputEvaluable, typename std::enable_if<!std::is_arithmetic<typename InputEvaluable::matrix_type>::value>::type>;

    template <class InputEvaluable, class OutputEvaluable>
    class BlockEvaluable<InputEvaluable, OutputEvaluable, typename std::enable_if<std::is_arithmetic<typename InputEvaluable::matrix_type>::value>::type>;

    template <class LeftEvaluable, class RightEvaluable>
    class CastEvaluable;

    template <typename EvaluableT>
    class SkewEvaluable;

    template <typename EvaluableT>
    class TransposeEvaluable;

    template <typename EvaluableT>
    class VeeEvaluable;

    template <typename EvaluableT, int rowsNumber = -1>
    class ConstructorByRows;

    template <typename EvaluableT, int colsNumber = -1>
    class ConstructorByCols;

    template <typename EvaluableT>
    class VariableFromExpressionEvaluable;

    template <typename CompositeEvaluable, typename leftEvaluable, typename rightEvaluable>
    class HorzcatEvaluable;

    template <typename CompositeEvaluable, typename TopEvaluable, typename BottomEvaluable>
    class VertcatEvaluable;

    template<typename EvaluableT>
    class SqueezeEvaluable;

    template<typename EvaluableT>
    class AutogeneratedEvaluable;

    template<typename EvaluableT>
    class TreeComponent;

    template<typename EvaluableT>
    class MultipleExpressionsEvaluator;

    template<typename EvaluableT>
    class MultipleCompiledExpressions;

    template<typename EvaluableT>
    class MultipleSqueezedExpressions;

    template<size_t startIndex, typename... Evaluables>
    class DependenciesSet;

    template<size_t startIndex, typename Evaluable>
    class DependenciesSet<startIndex, Evaluable>;

    template<typename EvaluableT>
    class DependenciesSet<0, EvaluableT>;

    template<size_t startIndex, class Evaluable0, class... OtherEvaluables>
    class DependenciesSet<startIndex, Evaluable0, OtherEvaluables...>;

    template<class Evaluable0, class... OtherEvaluables>
    class DependenciesSet<0, Evaluable0, OtherEvaluables...>;

    template<>
    class DependenciesSet<0>;

    template<typename... Evaluables>
    DependenciesSet<0, Evaluables...> make_dependencies_set(ExpressionComponent<Evaluables>&& ...dependencies);

    template<typename... Evaluables>
    DependenciesSet<0, Evaluables...> make_dependencies_set(const ExpressionComponent<Evaluables>& ...dependencies);

    template<typename EvaluableT>
    class AddendsExpander;

    template<typename GenericsMatrix, typename OutputMatrix>
    class CompiledEvaluable;


    /** Useful typedefs
     **/

    typedef Evaluable<LEVI_DEFAULT_MATRIX_TYPE> DefaultEvaluable;

    template<int rows, int cols>
    using DefaultFixedSizeEvaluable = Evaluable<LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, cols)>;

    typedef ExpressionComponent<DefaultEvaluable> Expression;

    typedef ExpressionComponent<Evaluable<double>> ScalarExpression;

    template<int rows, int cols>
    using FixedSizeExpression = ExpressionComponent<DefaultFixedSizeEvaluable<rows, cols>>;

    typedef ExpressionComponent<Evaluable<LEVI_DEFAULT_VECTOR_TYPE>> ColumnExpression;

    template<int rows, int cols>
    using FixedSizeColumnExpression = ExpressionComponent<Evaluable<LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, 1)>>;

    typedef EvaluableVariable<LEVI_DEFAULT_VECTOR_TYPE> DefaultVariableEvaluable;

    template <int dimension>
    using DefaultFixedSizeVariableEvaluable = EvaluableVariable<LEVI_DEFAULT_MATRIX_FIX_TYPE(dimension, 1)>;

    typedef ExpressionComponent<DefaultVariableEvaluable> Variable;

    template <int dimension>
    using FixedSizeVariable = ExpressionComponent<DefaultFixedSizeVariableEvaluable<dimension>>;

    typedef ExpressionComponent<EvaluableVariable<double>> ScalarVariable;

    typedef ExpressionComponent<ConstantEvaluable<LEVI_DEFAULT_MATRIX_TYPE>> Constant;

    typedef ExpressionComponent<MutableEvaluable<LEVI_DEFAULT_MATRIX_TYPE>> Mutable;

    template<int rows, int cols>
    using FixedSizeConstant = ExpressionComponent<ConstantEvaluable<LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, cols)>>;

    typedef ExpressionComponent<IdentityEvaluable<LEVI_DEFAULT_MATRIX_TYPE>> Identity;

    template<int rows, int cols>
    using FixedSizeIdentity = ExpressionComponent<IdentityEvaluable<LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, cols)>>;

    typedef ExpressionComponent<NullEvaluable<LEVI_DEFAULT_MATRIX_TYPE>> Null;

    template<int rows, int cols>
    using FixedSizeNull = ExpressionComponent<NullEvaluable<LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, cols)>>;

    typedef ExpressionComponent<ConstantEvaluable<double>> Scalar;

    typedef ExpressionComponent<MutableEvaluable<double>> ScalarMutable;

}

#endif // LEVI_FORWARDDECLARATIONS_H
