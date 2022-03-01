/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_VARIABLE_H
#define LEVI_VARIABLE_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/VariableBase.h>
#include <levi/BasicEvaluables.h>

/**
 * @brief The EvaluableVariable class, specialized for vector values.
 *
 * This class inherits from both VariableBase and Evaluable. This allows to use the variables inside expressions as they were evaluables,
 * thus allowing for seamless application of arithmetic operators.
 */
template <typename Vector>
class levi::EvaluableVariable<Vector, typename std::enable_if<!std::is_arithmetic<Vector>::value>::type> : public levi::VariableBase,
                                                                                                           public levi::Evaluable<Vector>
{
protected:
    levi::ExpressionComponent<levi::IdentityEvaluable<typename levi::Evaluable<Vector>::derivative_evaluable::matrix_type>> m_identityDerivative;

    template<typename OtherVector, bool isScalar>
    void copy_constant(bool_value<isScalar>, const OtherVector& rhs);

    template<typename OtherVector>
    void copy_constant(bool_value<true>, const OtherVector& rhs) {
        static_assert (this->rows_at_compile_time == Eigen::Dynamic || this->rows_at_compile_time == 1, "Cannot copy a scalar to this variable.");

        this->m_evaluationBuffer.resize(1,1);

        this->m_evaluationBuffer(0,0) = rhs;
    }

    template<typename OtherVector>
    void copy_constant(bool_value<false>, const OtherVector& rhs) {
        static_assert (OtherVector::ColsAtCompileTime == 1, "The chosen VectorType for the rhs should have exactly one column at compile time.");
        assert(rhs.size() == this->dimension());
        this->m_evaluationBuffer = rhs;
    }

public:

    EvaluableVariable(Eigen::Index dimension, const std::string& name)
        : levi::VariableBase(dimension, name)
        , levi::Evaluable<Vector>(dimension, 1, name)
        , m_identityDerivative(dimension, dimension, "d " + name + "/(d " + name + ")")
    {
        static_assert (Vector::ColsAtCompileTime == 1, "The chosen VectorType for the Variable should have exactly one column at compile time.");

        this->m_evaluationBuffer.setZero();

    }

    template <typename otherVector>
    EvaluableVariable(const EvaluableVariable<otherVector>& other) = delete;

    template <typename otherVector>
    EvaluableVariable(EvaluableVariable<otherVector>&& other) = delete;

    virtual ~EvaluableVariable() override;

    /**
     * @brief Assignement operator to set the values of the variable.
     *
     * The right hand side has to be a vector.
     */
    template<typename otherVector>
    void operator=(const otherVector& rhs) {
        copy_constant(bool_value<std::is_arithmetic<otherVector>::value>(), rhs);
        this->resetEvaluationRegister();
    }

    virtual const Vector& evaluate() override {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Vector>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                     std::shared_ptr<levi::VariableBase> variable) override {
        levi::unused(column);
        assert(column == 0);
        if ((this->variableName() == variable->variableName()) && (this->dimension() == variable->dimension())) {
            return m_identityDerivative;
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Vector>::derivative_evaluable::matrix_type>>(this->dimension(), variable->dimension(),
                                                                                                                                       "d " + variableName() + "/(d " + variable->variableName() + ")");
        }
    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) override {
        return ((this->variableName() == variable->variableName()) && (this->dimension() == variable->dimension()));
    }

    virtual std::vector<levi::Registrar*> getDependencies() override {
        return {this};
    }

};
template <typename Vector>
levi::EvaluableVariable<Vector, typename std::enable_if<!std::is_arithmetic<Vector>::value>::type>::~EvaluableVariable() { }

/**
 * @brief The EvaluableVariable class, specialized for scalar values.
 *
 * This class inherits from both VariableBase and Evaluable. This allows to use the variables inside expressions as they were evaluables,
 * thus allowing for seamless application of arithmetic operators.
 */
template <typename Scalar>
class levi::EvaluableVariable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type> : public levi::VariableBase,
                                                                                                          public levi::Evaluable<Scalar>
{
protected:
    levi::ExpressionComponent<levi::IdentityEvaluable<typename levi::Evaluable<Scalar>::derivative_evaluable::matrix_type>> m_identityDerivative;

public:

    EvaluableVariable(Eigen::Index dimension, const std::string& name)
        : levi::VariableBase(1, name)
          , levi::Evaluable<Scalar>(1, 1, name)
          , m_identityDerivative(1, 1, "d " + name + "/(d " + name + ")")
    {
        assert(dimension == 1 && "double has been choosen as a type but a dimension different from 1 has been used.");
    }

    EvaluableVariable(const std::string& name)
        : levi::VariableBase(1, name)
        , levi::Evaluable<Scalar>(0, name)
        , m_identityDerivative(1,1)
    { }

    template <typename otherVector>
    EvaluableVariable(const EvaluableVariable<otherVector>& other) = delete;

    template <typename otherVector>
    EvaluableVariable(EvaluableVariable<otherVector>&& other) = delete;

    virtual ~EvaluableVariable() override;

    /**
     * @brief Assignement operator to set the values of the variable.
     *
     * The right hand side has to be a scalar.
     */
    void operator=(const Scalar& rhs) {
        this->m_evaluationBuffer = rhs;
        this->resetEvaluationRegister();
    }

    virtual const Scalar& evaluate() override {
        return this->m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::Evaluable<Scalar>::derivative_evaluable> getNewColumnDerivative(Eigen::Index column,
                                                                                                                     std::shared_ptr<levi::VariableBase> variable) override {
        levi::unused(column);
        assert(column == 0);
        if ((this->variableName() == variable->variableName()) && (this->dimension() == variable->dimension())) {
            return m_identityDerivative;
        } else {
            return levi::ExpressionComponent<levi::NullEvaluable<typename levi::Evaluable<Scalar>::derivative_evaluable::matrix_type>>(this->dimension(), variable->dimension());
        }
    }

    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) override {
        return ((this->variableName() == variable->variableName()) && (this->dimension() == variable->dimension()));
    }

    virtual std::vector<levi::Registrar*> getDependencies() override {
        return {this};
    }

};
template <typename Scalar>
levi::EvaluableVariable<Scalar, typename std::enable_if<std::is_arithmetic<Scalar>::value>::type>::~EvaluableVariable() { }

#endif // LEVI_VARIABLE_H
