/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_TREEEXPANDER_H
#define LEVI_TREEEXPANDER_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Expression.h>
#include <levi/TypeDetector.h>
#include <map>

namespace levi {

    typedef std::map<levi::EvaluableType, std::vector<size_t>> AddedExpressions;

    template<typename EvaluableT>
    size_t expandTree(const levi::ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic>>>& node,
                      bool expandForSqueeze,
                      std::vector<levi::TreeComponent<EvaluableT>>& expandedExpression,
                      std::vector<size_t>& generics);

    template<typename EvaluableT>
    size_t expandTree(const levi::ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic>>>& node,
                      bool expandForSqueeze,
                      std::vector<levi::TreeComponent<EvaluableT>>& expandedExpression,
                      std::vector<size_t>& generics, AddedExpressions& alreadyAdded);
}

template<typename EvaluableT>
class levi::TreeComponent {

    using Type = levi::EvaluableType;

    Eigen::Index m_rows;
    Eigen::Index m_cols;

public:

    typedef Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic> SqueezedMatrix;

    levi::ExpressionComponent<levi::Evaluable<SqueezedMatrix>> partialExpression;

    Type type;

    SqueezedMatrix buffer;

    size_t lhsIndex;
    size_t rhsIndex;

    levi::BlockType block;
    typename EvaluableT::value_type exponent;


    TreeComponent(const levi::ExpressionComponent<levi::Evaluable<SqueezedMatrix>>& expression, bool expandForSqueeze)
        : partialExpression(expression)
          , type(expression.info().type)
          , m_rows(expression.rows())
          , m_cols(expression.cols())
          , lhsIndex(0)
          , rhsIndex(0)
    {
        if (!expandForSqueeze && (type == Type::Null || type == Type::Identity || type == Type::Constant || type == Type::Horzcat || type == Type::Vertcat)) {
            type = Type::Generic; //This is supposed to be a temporary fix for evaluables not yet supported by autogeneration
        }

        if (expandForSqueeze || type == Type::Generic) {
            buffer.resize(expression.rows(), expression.cols());

            if (type == Type::Null) {
                buffer.setZero();
            }

            if (type == Type::Identity) {
                buffer.setIdentity();
            }
        }

        if (type == Type::Constant) {
            buffer = partialExpression.evaluate();
        }

        if (type != Type::Generic) {
            block = expression.info().block;
            exponent = expression.info().exponent;
        }
    }

    Eigen::Index rows() const {
        return m_rows;
    }

    Eigen::Index cols() const {
        return m_cols;
    }

};

template<typename EvaluableT>
size_t levi::expandTree(const levi::ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic>>>& node,
                        bool expandForSqueeze,
                        std::vector<levi::TreeComponent<EvaluableT>>& expandedExpression,
                        std::vector<size_t>& generics) {
    levi::AddedExpressions alreadyAdded;
    return levi::expandTree(node, expandForSqueeze, expandedExpression, generics, alreadyAdded);
}

template<typename EvaluableT>
size_t levi::expandTree(const levi::ExpressionComponent<levi::Evaluable<Eigen::Matrix<typename EvaluableT::value_type, Eigen::Dynamic, Eigen::Dynamic>>>& node,
                        bool expandForSqueeze,
                        std::vector<levi::TreeComponent<EvaluableT>>& expandedExpression,
                        std::vector<size_t>& generics, levi::AddedExpressions& alreadyAdded) {

    using Type = levi::EvaluableType;

    //Returns the position in which the node has been saved
    assert(node.isValidExpression());

    levi::EvaluableType type;

    levi::TreeComponent<EvaluableT> newComponent(node, expandForSqueeze);

    type = newComponent.type;

    size_t currentIndex = expandedExpression.size();

    if (type == Type::Generic) {
        for (size_t generic : generics) {
            if (expandedExpression[generic].partialExpression == newComponent.partialExpression) {
                // This is the case where the same Generic has been added before.
                return generic;
            }
        }
        expandedExpression.push_back(newComponent);
        currentIndex = expandedExpression.size() - 1;
        generics.push_back(currentIndex);
        return currentIndex;
    }

    if (type == Type::Null || type == Type::Identity || type == Type::Constant) {
        expandedExpression.push_back(newComponent);
        currentIndex = expandedExpression.size() - 1;
        return currentIndex;
    }

    if (type == Type::Sum || type == Type::Subtraction || type == Type::Product ||
        type == Type::Division || type == Type::Vertcat || type == Type::Horzcat) {
        size_t lhsIndex = expandTree(node.info().lhs, expandForSqueeze, expandedExpression, generics, alreadyAdded);
        newComponent.lhsIndex = lhsIndex;
        size_t rhsIndex = expandTree(node.info().rhs, expandForSqueeze, expandedExpression, generics, alreadyAdded);
        newComponent.rhsIndex = rhsIndex;

        if (lhsIndex < currentIndex && rhsIndex < currentIndex) { //This means that both the lhs and the rhs were added before
            bool isNew = true;
            size_t i = 0;

            std::vector<size_t>& others = alreadyAdded[type];

            while (isNew && i < others.size()) {
                isNew = isNew && (expandedExpression[others[i]].partialExpression != newComponent.partialExpression);
                ++i;
            }

            if (!isNew) {
                return others[i-1];
            }

        }

        expandedExpression.push_back(newComponent);
        currentIndex = expandedExpression.size() - 1;
        alreadyAdded[type].push_back(currentIndex);
        return currentIndex;
    }

    if (type == Type::InvertedSign || type == Type::Pow || type == Type::Transpose || type == Type::Row ||
        type == Type::Column || type == Type::Element || type == Type::Block) {
        size_t lhsIndex = expandTree(node.info().lhs, expandForSqueeze, expandedExpression, generics, alreadyAdded);
        newComponent.lhsIndex = lhsIndex;

        if (lhsIndex < currentIndex) { //This means that the lhs was already before
            bool isNew = true;
            size_t i = 0;

            std::vector<size_t>& others = alreadyAdded[type];

            while (isNew && i < others.size()) {
                isNew = isNew && (expandedExpression[others[i]].partialExpression != newComponent.partialExpression);
                ++i;
            }

            if (!isNew) {
                return others[i-1];
            }

        }

        expandedExpression.push_back(newComponent);
        currentIndex = expandedExpression.size() - 1;
        alreadyAdded[type].push_back(currentIndex);
        return currentIndex;
    }

    assert(false && "Case not considered.");
    return expandedExpression.size();
}

#endif // LEVI_TREEEXPANDER_H
