/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_TYPEDETECTOR_H
#define LEVI_TYPEDETECTOR_H

#include <levi/ForwardDeclarations.h>

namespace levi {
    enum class EvaluableType {
        Generic,
        Null,
        Identity,
        Constant,
        Sum,
        Subtraction,
        InvertedSign,
        Product,
        Division,
        Pow,
        Transpose,
        Row,
        Column,
        Element,
        Block,
        Vertcat,
        Horzcat
    };

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<EvaluableT>&) {
        return EvaluableType::Generic;
    }

    template <typename Matrix>
    EvaluableType detectType(const ExpressionComponent<NullEvaluable<Matrix>>&) {
        return EvaluableType::Null;
    }

    template <typename Matrix>
    EvaluableType detectType(const ExpressionComponent<IdentityEvaluable<Matrix>>&) {
        return EvaluableType::Identity;
    }

    template <typename EvaluableLhs, typename EvaluableRhs>
    EvaluableType detectType(const ExpressionComponent<SumEvaluable<EvaluableLhs, EvaluableRhs>>&) {
        return EvaluableType::Sum;
    }

    template <typename EvaluableLhs, typename EvaluableRhs>
    EvaluableType detectType(const ExpressionComponent<SubtractionEvaluable<EvaluableLhs, EvaluableRhs>>&) {
        return EvaluableType::Subtraction;
    }

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<SignInvertedEvaluable<EvaluableT>>&) {
        return EvaluableType::InvertedSign;
    }

    template <typename EvaluableLhs, typename EvaluableRhs>
    EvaluableType detectType(const ExpressionComponent<ProductEvaluable<EvaluableLhs, EvaluableRhs>>&) {
        return EvaluableType::Product;
    }

    template <typename EvaluableLhs, typename EvaluableRhs>
    EvaluableType detectType(const ExpressionComponent<DivisionEvaluable<EvaluableLhs, EvaluableRhs>>&) {
        return EvaluableType::Division;
    }

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<PowEvaluable<EvaluableT>>&) {
        return EvaluableType::Pow;
    }

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<TransposeEvaluable<EvaluableT>>&) {
        return EvaluableType::Transpose;
    }

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<RowEvaluable<EvaluableT>>&) {
        return EvaluableType::Row;
    }

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<ColEvaluable<EvaluableT>>&) {
        return EvaluableType::Column;
    }

    template <typename EvaluableT>
    EvaluableType detectType(const ExpressionComponent<ElementEvaluable<EvaluableT>>&) {
        return EvaluableType::Element;
    }

    template <typename EvaluableT, typename EvaluableOut>
    EvaluableType detectType(const ExpressionComponent<BlockEvaluable<EvaluableT, EvaluableOut>>&) {
        return EvaluableType::Block;
    }

}

#endif // LEVI_TYPEDETECTOR_H
