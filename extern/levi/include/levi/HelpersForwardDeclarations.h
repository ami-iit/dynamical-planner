/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_HELPERSFORWARDDECLARATIONS_H
#define LEVI_HELPERSFORWARDDECLARATIONS_H

#include <Eigen/Core>
#include <vector>
#include <memory>
#include <string>
#include <cassert>
#include <cmath>
#include <unordered_map>
#include <mutex>

#ifndef LEVI_DEFAULT_MATRIX_TYPE
#define LEVI_DEFAULT_MATRIX_TYPE Eigen::MatrixXd
#endif

#ifndef LEVI_DEFAULT_MATRIX_FIX_TYPE
#define LEVI_DEFAULT_MATRIX_FIX_TYPE(rows, cols) Eigen::Matrix<double, rows, cols>
#endif

#ifndef LEVI_DEFAULT_VECTOR_TYPE
#define LEVI_DEFAULT_VECTOR_TYPE Eigen::VectorXd
#endif

namespace levi {

    template <typename... Args> inline void unused(Args&&...) {}

    template<bool T>
    struct bool_value { };

    template<size_t value>
    struct size_t_value { };

    template <typename Scalar_lhs, typename Scalar_rhs>
    struct scalar_sum_return;

    template<int lhsRows, int lhsCols, int rhsRows, int rhsCols, class Enabler = void>
    struct is_valid_sum : std::false_type {};

    template<int lhsRows, int lhsCols, int rhsRows, int rhsCols>
    struct is_valid_sum<lhsRows, lhsCols, rhsRows, rhsCols,
            typename std::enable_if<((lhsRows == rhsRows) || (lhsRows == Eigen::Dynamic) || (rhsRows == Eigen::Dynamic)) &&
    ((lhsCols == rhsCols) || (lhsCols == Eigen::Dynamic) || (rhsCols == Eigen::Dynamic))>::type> : std::true_type {};

    template <typename Matrix_lhs, typename Matrix_rhs, class Enabler = void>
    struct matrix_sum_return;

    template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_sum_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
            typename std::enable_if<is_valid_sum<lhsRows, lhsCols, rhsRows, rhsCols>::value>::type>;

    template<typename Scalar_lhs, typename Scalar_rhs>
    struct matrix_sum_return<Scalar_lhs, Scalar_rhs,
            typename std::enable_if<std::is_arithmetic<Scalar_lhs>::value && std::is_arithmetic<Scalar_rhs>::value>::type>;

    template<typename Scalar, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_sum_return<Scalar, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
            typename std::enable_if<std::is_arithmetic<Scalar>::value && is_valid_sum<1,1, rhsRows, rhsCols>::value>::type>;

    template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar>
    struct matrix_sum_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Scalar,
            typename std::enable_if<std::is_arithmetic<Scalar>::value && is_valid_sum<1,1, lhsRows, lhsCols>::value>::type>;

    template<int lhsRows, int lhsCols, int rhsRows, int rhsCols, class Enabler = void>
    struct is_valid_product : std::false_type {};

    template<int lhsRows, int lhsCols, int rhsRows, int rhsCols>
    struct is_valid_product<lhsRows, lhsCols, rhsRows, rhsCols,
                            typename std::enable_if<((lhsRows == 1 || lhsRows == Eigen::Dynamic) && (lhsCols == 1 || lhsCols == Eigen::Dynamic)) ||
                                                    ((rhsRows == 1 || rhsRows == Eigen::Dynamic) && (rhsCols == 1 || rhsCols == Eigen::Dynamic)) ||
                                                    lhsCols == Eigen::Dynamic || rhsRows == Eigen::Dynamic || lhsCols == rhsRows>::type> : std::true_type {};

    template <typename Scalar_lhs, typename Scalar_rhs>
    struct scalar_product_return;

    template <typename Matrix_lhs, typename Matrix_rhs, class Enabler = void>
    struct matrix_product_return;

    template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
            typename std::enable_if<is_valid_product<lhsRows, lhsCols, rhsRows, rhsCols>::value &&
            !(lhsRows == 1 && lhsCols == 1) && !(rhsRows == 1 && rhsCols == 1)>::type>;

    template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
                                 typename std::enable_if<lhsRows == 1 && lhsCols == 1 && (rhsRows != 1 || rhsCols != 1)>::type>;

    template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
            typename std::enable_if<rhsRows == 1 && rhsCols == 1 && (lhsRows != 1 || lhsCols != 1)>::type>;

    template<typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
                                 typename std::enable_if<rhsRows == 1 && rhsCols == 1 && lhsRows == 1 && lhsCols == 1>::type>;

    template<typename Scalar, typename Scalar_rhs, int rhsRows, int rhsCols, int rhsOptions, int rhsMaxRows, int rhsMaxCols>
    struct matrix_product_return<Scalar, Eigen::Matrix<Scalar_rhs, rhsRows, rhsCols, rhsOptions, rhsMaxRows, rhsMaxCols>,
            typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template<typename Scalar, typename Scalar_lhs, int lhsRows, int lhsCols, int lhsOptions, int lhsMaxRows, int lhsMaxCols>
    struct matrix_product_return<Eigen::Matrix<Scalar_lhs, lhsRows, lhsCols, lhsOptions, lhsMaxRows, lhsMaxCols>, Scalar,
            typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template<typename Scalar_lhs, typename Scalar_rhs>
    struct matrix_product_return<Scalar_lhs, Scalar_rhs,
            typename std::enable_if<std::is_arithmetic<Scalar_lhs>::value && std::is_arithmetic<Scalar_rhs>::value>::type>;

    template<typename Matrix, class Enabler = void>
    struct dynamic_block_return;

    template<typename Matrix>
    struct dynamic_block_return<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template<typename Scalar>
    struct dynamic_block_return<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template<typename Matrix, int rows, int cols, class Enabler = void>
    struct fixed_block_return;

    template<typename Matrix, int rows, int cols>
    struct fixed_block_return<Matrix, rows, cols, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>;

    template<typename Scalar, int rows, int cols>
    struct fixed_block_return<Scalar, rows, cols, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>;

    template<typename EvaluableT, class Enabler = void>
    struct transpose_type;

    template<typename EvaluableT>
    struct transpose_type<EvaluableT, typename std::enable_if<!std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template<typename EvaluableT>
    struct transpose_type<EvaluableT, typename std::enable_if<std::is_arithmetic<typename EvaluableT::matrix_type>::value>::type>;

    template<typename Scalar>
    struct TripletStruct;

    struct BlockType {
        Eigen::Index rows = -1;
        Eigen::Index cols = -1;
        Eigen::Index startRow = -1;
        Eigen::Index startCol = -1;

        bool operator==(const BlockType& other) const {
            return (rows == other.rows) && (cols == other.cols) && (startRow == other.startRow) && (startCol == other.startCol);
        }
    };

    template<typename Scalar>
    using Triplet = TripletStruct<Scalar>;

}

#endif // LEVI_HELPERSFORWARDDECLARATIONS_H
