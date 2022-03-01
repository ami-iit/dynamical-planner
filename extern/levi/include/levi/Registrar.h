/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_REGISTRAR_H
#define LEVI_REGISTRAR_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Expression.h>

class levi::Registrar {
protected:

    /**
     * @brief Register to keep track of which parent read the new values.
     */
    std::vector<int> m_evaluationRegister;

    int m_isNewCounter;

    /**
     * @brief Register to keep track of which index can be reused
     */
    std::vector<size_t> m_IDRecycleBin;

    /**
     * @brief Reset the evaluation register to notify of new values
     */
    void resetEvaluationRegister() {
        m_isNewCounter++;
        m_alreadyComputed = false;
    }

    /**
     * @brief Check if the buffer has been filled at least once since the evaluation register has been resetted
     */
    bool m_alreadyComputed;

    struct ValidityChecker { };

    std::shared_ptr<ValidityChecker> m_thisIsValid;

    class Dependency {
        Registrar* m_ptr;
        size_t m_id;
        std::weak_ptr<ValidityChecker> m_isValid;

    public:

        Dependency() = delete;

        Dependency(const Dependency& other) {
            if (isAlive()) {
                m_ptr->deleteID(m_id);
            }

            m_ptr = other.m_ptr;
            assert(m_ptr);
            m_isValid = m_ptr->getValidityChecker();
            assert(isAlive());
            m_id = m_ptr->getNewCallerID();
        }

        Dependency(Dependency&& other) {
            if (isAlive()) {
                m_ptr->deleteID(m_id);
            }

            m_ptr = other.m_ptr;
            assert(m_ptr);
            m_isValid = m_ptr->getValidityChecker();
            assert(isAlive());
            m_id = m_ptr->getNewCallerID();
        }

        Dependency(Registrar* other)
        {
            m_ptr = other;
            assert(other);
            m_isValid = other->getValidityChecker();
            assert(isAlive());
            m_id = other->getNewCallerID();
        }

        ~Dependency() {
            if (isAlive()) {
                m_ptr->deleteID(m_id);
            }
        }

        bool operator!=(Registrar* other) const {
            assert(isAlive() && "The dependency was deleted while it was still needed.");
            return (m_ptr != other);
        }

        Registrar* get() {
            assert(isAlive() && "The dependency was deleted while it was still needed.");
            return m_ptr;
        }

        bool isNew() const {
            assert(isAlive() && "The dependency was deleted while it was still needed.");
            return m_ptr->isNew(m_id);
        }

        bool isAlive() const {
            return !m_isValid.expired();
        }

        void check() {
            assert(isAlive() && "The dependency was deleted while it was still needed.");
            m_ptr->checkID(m_id);
        }
    };

    /**
     * @brief Vector storing pointers to the dependencies registrar, plus the EvaluationID
     *
     * @warning Here we assume that the ptr will be always be alive when evaluating. Infact,
     * if it is a dependecy, I can expect the evaluable to save a pointer to the expression. Indeed, dependencies
     * cannot be evaluated via the pointer to the registrar, thus another pointer is needed to exploit their values.
     */
    std::vector<Dependency> m_dependencies;

    void addDependencies(const std::vector<Registrar*>& dependencies) {
        for (auto& dep : dependencies) {
            assert(dep);
            bool isNewDep = true;
            size_t i = 0;

            while (isNewDep && i < m_dependencies.size()) {
                isNewDep = isNewDep && m_dependencies[i] != dep;
                ++i;
            }

            if (isNewDep) {
                m_dependencies.emplace_back(dep);
            }
        }
    }

    template<typename... OtherVectors>
    void addDependencies(const std::vector<Registrar*>& firstDependencies, const OtherVectors& ...otherDependencies) {
        addDependencies(firstDependencies);
        addDependencies(otherDependencies...);
    }

    template <typename EvaluableT>
    void addDependencies(const levi::ExpressionComponent<EvaluableT>& dependentExpression) {
        for (auto& dep : dependentExpression.getDependencies()) {
            assert(dep);
            bool isNewDep = true;
            size_t i = 0;

            while (isNewDep && i < m_dependencies.size()) {
                isNewDep = isNewDep && m_dependencies[i] != dep;
                ++i;
            }

            if (isNewDep) {
                m_dependencies.emplace_back(dep);
            }
        }
    }

    template<typename EvaluableT, typename... OtherExpressions>
    void addDependencies(const levi::ExpressionComponent<EvaluableT>& firstDependencies, const OtherExpressions& ...otherDependencies) {
        addDependencies(firstDependencies);
        addDependencies(otherDependencies...);
    }

public:

    Registrar()
        : m_isNewCounter(0)
          , m_alreadyComputed(false)
          , m_thisIsValid(std::make_shared<ValidityChecker>())
    { }

    virtual ~Registrar() { }

    /**
     * @brief Check if the evaluable has a new value which has not been read by the caller
     * @param callerID The ID of caller
     * @return True if new
     */
    bool isNew(size_t callerID) {

        bool isNew = false;
        size_t i = 0;

        while (!isNew && i < m_dependencies.size()) {
            isNew = isNew || (m_dependencies[i] != this && m_dependencies[i].isNew());
            ++i;
        }

        if (isNew) {
            resetEvaluationRegister();
        }

        return m_evaluationRegister[callerID] != m_isNewCounter;
    }

    /**
     * @brief Check whether the evaluable depends on a specified variable
     * @param variable The variable of interest
     * @return True if dependent
     *
     */
    virtual bool isDependentFrom(std::shared_ptr<levi::VariableBase> variable) {
        bool isDependent = false;
        size_t i = 0;
        while (!isDependent && i < m_dependencies.size()) {
            isDependent = isDependent || (m_dependencies[i] != this && m_dependencies[i].get()->isDependentFrom(variable));
            ++i;
        }

        return isDependent;
    }

    virtual std::vector<levi::Registrar*> getDependencies() {
        std::vector<Registrar*> deps;
        for (auto& dep : m_dependencies) {
            deps.push_back(dep.get());
        }
        return deps;
    }

    void dependenciesChecked() {
        for (auto& dep : m_dependencies) {
            dep.check();
        }
    }

    void checkID(size_t callerID) {
        m_evaluationRegister[callerID] = m_isNewCounter;
    }

    /**
     * @brief Get a new ID to check whether the evaluable has a new value
     * @return An ID to use to check whether the evaluable has already been called since new data has come in.
     */
    size_t getNewCallerID() {
        if (m_IDRecycleBin.size()) {
            size_t newIndex = m_IDRecycleBin.back();
            m_evaluationRegister[newIndex] = m_isNewCounter - 1;
            m_IDRecycleBin.pop_back();
            return newIndex;
        }
        m_evaluationRegister.push_back(m_isNewCounter - 1);
        return m_evaluationRegister.size() -1;
    }

    std::weak_ptr<ValidityChecker> getValidityChecker() const {
        return m_thisIsValid;
    }


    /**
     * @brief Delete a previousy created ID
     * @param index The index to remove
     * The index will be contained in the buffer m_IDRecycleBin to reuse previously deleted IDs. This avoids the
     * m_evaluationRegister to grow indefinitely
     */
    void deleteID(size_t index) {
        m_evaluationRegister[index] = m_isNewCounter - 1;
        m_IDRecycleBin.push_back(index);
    }

};

#endif // LEVI_REGISTRAR_H
