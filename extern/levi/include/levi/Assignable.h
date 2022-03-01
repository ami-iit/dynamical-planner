/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_ASSIGNABLE_H
#define LEVI_ASSIGNABLE_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Evaluable.h>
#include <levi/VariableBase.h>

/**
 * @brief The Assignable evaluable
 *
 * Evaluable which can be assigned through the operator =
 * Specialization for matrix types
 *
 */
template <typename Matrix>
class levi::Assignable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type> : public levi::Evaluable<Matrix>
{
public:

    Assignable(std::string name)
        : levi::Evaluable<Matrix>(name)
    {
        this->addDependencies({this});
    }

    Assignable(const Matrix& initialValue, std::string name)
        : levi::Evaluable<Matrix>(initialValue, name)
    {
        this->addDependencies({this});
    }

    Assignable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Evaluable<Matrix>(rows, cols, name)
    {
        this->addDependencies({this});
    }

    virtual ~Assignable() override;

    template <typename MatrixRhs>
    void operator=(const MatrixRhs& rhs) {
        this->m_evaluationBuffer = rhs;
        this->resetEvaluationRegister();
    }

};
template <typename Matrix>
levi::Assignable<Matrix, typename std::enable_if<!std::is_arithmetic<Matrix>::value>::type>::~Assignable() { }

/**
 * @brief The Assignable evaluable
 *
 * Evaluable which can be assigned through the operator =
 * Specialization for scalar types
 *
 */
template <typename Scalar>
class levi::Assignable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::Evaluable<Scalar>
{
public:

    Assignable(const std::string& name)
        : levi::Evaluable<Scalar>(name)
    {
        this->addDependencies(this);
    }

    Assignable(Eigen::Index rows, Eigen::Index cols, const std::string& name)
        : levi::Evaluable<Scalar>(rows, cols, name)
    {
        this->addDependencies(this);
    }


    Assignable(const Scalar& initialValue)
        : levi::Evaluable<Scalar>(initialValue)
    {
        this->addDependencies(this);
    }

    Assignable(const Scalar& initialValue, const std::string& name)
        : levi::Evaluable<Scalar>(initialValue, name)
    {
        this->addDependencies(this);
    }

    virtual ~Assignable() override;

    void operator=(const Scalar& rhs) {
        this->m_evaluationBuffer = rhs;
        this->resetEvaluationRegister();
    }

};
template <typename Scalar>
levi::Assignable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~Assignable() { }

#endif // LEVI_ASSIGNABLE_H
