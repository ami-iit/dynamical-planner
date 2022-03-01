/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef LEVI_DEPENDENCIESSET_H
#define LEVI_DEPENDENCIESSET_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>
#include <levi/Expression.h>

template<size_t startIndex, typename EvaluableT>
class levi::DependenciesSet<startIndex, EvaluableT> {
    levi::ExpressionComponent<EvaluableT> m_dependency;

public:

    DependenciesSet(const levi::ExpressionComponent<EvaluableT>& dependency)
        : m_dependency(dependency)
    { }

    DependenciesSet(levi::ExpressionComponent<EvaluableT>&& dependency)
        : m_dependency(std::forward<levi::ExpressionComponent<EvaluableT>>(dependency))
    { }

    DependenciesSet(const DependenciesSet<startIndex, EvaluableT>& other)
        : m_dependency(other.m_dependency)
    { }

    template<size_t index>
    const levi::ExpressionComponent<EvaluableT>& get(levi::size_t_value<index>) const;

    const levi::ExpressionComponent<EvaluableT>& get(levi::size_t_value<startIndex>) const {
        return m_dependency;
    }

    template<size_t index>
    bool isNew(levi::size_t_value<index>) const;

    bool isNew(levi::size_t_value<startIndex>) const {
        bool isNew = m_dependency.isNew();
        return isNew;
    }

    template<size_t index>
    const typename EvaluableT::matrix_type& evaluate(levi::size_t_value<index>);

    const typename EvaluableT::matrix_type& evaluate(levi::size_t_value<startIndex>) {
        return m_dependency.evaluate();
    }
};

template<typename EvaluableT>
class levi::DependenciesSet<0, EvaluableT> {
    levi::ExpressionComponent<EvaluableT> m_dependency;

public:

    DependenciesSet(const levi::ExpressionComponent<EvaluableT>& dependency)
        : m_dependency(dependency)
    { }

    DependenciesSet(levi::ExpressionComponent<EvaluableT>&& dependency)
        : m_dependency(dependency)
    { }

    DependenciesSet(const DependenciesSet<0, EvaluableT>& other)
        : m_dependency(other.m_dependency)
    { }

    template<size_t index>
    const levi::ExpressionComponent<EvaluableT>& get(levi::size_t_value<index>) const;

    const levi::ExpressionComponent<EvaluableT>& get(levi::size_t_value<0>) const {
        return m_dependency;
    }

    template<size_t index>
    bool isNew(levi::size_t_value<index>) const;

    bool isNew(levi::size_t_value<0>) const {
        bool isNew = m_dependency.isNew();
        return isNew;
    }

    template<size_t index>
    const typename EvaluableT::matrix_type& evaluate(levi::size_t_value<index>);

    const typename EvaluableT::matrix_type& evaluate(levi::size_t_value<0>) {
        return m_dependency.evaluate();
    }
};

template<size_t startIndex, class Evaluable0, class... OtherEvaluables>
class levi::DependenciesSet<startIndex, Evaluable0, OtherEvaluables...> : public levi::DependenciesSet<startIndex, Evaluable0>,
                                                                          public levi::DependenciesSet<startIndex+1, OtherEvaluables...>
{
public:

    DependenciesSet(const levi::ExpressionComponent<Evaluable0>& firstDependency,
                    const levi::ExpressionComponent<OtherEvaluables>& ...otherDependencies)
        :levi::DependenciesSet<startIndex, Evaluable0>(firstDependency)
          ,levi::DependenciesSet<startIndex+1, OtherEvaluables...>(otherDependencies...)
    { }

    DependenciesSet(levi::ExpressionComponent<Evaluable0>&& firstDependency, levi::ExpressionComponent<OtherEvaluables>&& ...otherDependencies)
        :levi::DependenciesSet<startIndex, Evaluable0>(std::forward(firstDependency))
          ,levi::DependenciesSet<startIndex+1, OtherEvaluables...>(std::forward<levi::ExpressionComponent<OtherEvaluables>>(otherDependencies)...)
    { }

    using levi::DependenciesSet<startIndex, Evaluable0>::isNew;
    using levi::DependenciesSet<startIndex+1, OtherEvaluables...>::isNew;
    using levi::DependenciesSet<startIndex, Evaluable0>::get;
    using levi::DependenciesSet<startIndex+1, OtherEvaluables...>::get;
    using levi::DependenciesSet<startIndex, Evaluable0>::evaluate;
    using levi::DependenciesSet<startIndex+1, OtherEvaluables...>::evaluate;
};

template<class Evaluable0, class... OtherEvaluables>
class levi::DependenciesSet<0, Evaluable0, OtherEvaluables...> : public levi::DependenciesSet<0, Evaluable0>,
                                                                 public levi::DependenciesSet<1, OtherEvaluables...>
{

public:

    using levi::DependenciesSet<0, Evaluable0>::isNew;
    using levi::DependenciesSet<1, OtherEvaluables...>::isNew;
    using levi::DependenciesSet<0, Evaluable0>::get;
    using levi::DependenciesSet<1, OtherEvaluables...>::get;
    using levi::DependenciesSet<0, Evaluable0>::evaluate;
    using levi::DependenciesSet<1, OtherEvaluables...>::evaluate;

private:
    template<std::size_t currentIndex>
    bool areNewHelper(levi::size_t_value<currentIndex>) const {
        bool isFirstNew = isNew(levi::size_t_value<currentIndex>());
        return areNewHelper(levi::size_t_value<currentIndex + 1>()) || isFirstNew;
    }

    bool areNewHelper(levi::size_t_value<sizeof... (OtherEvaluables)>) const {
        return isNew(levi::size_t_value<sizeof... (OtherEvaluables)>());
    }

public:

    DependenciesSet(const levi::ExpressionComponent<Evaluable0>& firstDependency,
                    const levi::ExpressionComponent<OtherEvaluables>& ...otherDependencies)
        :levi::DependenciesSet<0, Evaluable0>(firstDependency)
          ,levi::DependenciesSet<1, OtherEvaluables...>(otherDependencies...)
    { }

    DependenciesSet(levi::ExpressionComponent<Evaluable0>&& firstDependency, levi::ExpressionComponent<OtherEvaluables>&& ...otherDependencies)
        :levi::DependenciesSet<0, Evaluable0>(std::forward(firstDependency))
          ,levi::DependenciesSet<1, OtherEvaluables...>(std::forward<levi::ExpressionComponent<OtherEvaluables>>(otherDependencies)...)
    { }

    size_t size() const {
        return sizeof... (OtherEvaluables) + 1;
    }

    bool areNew() const {
        return areNewHelper(levi::size_t_value<0>());
    }

};

template<>
class levi::DependenciesSet<0>
{
public:

    DependenciesSet()
    { }

    size_t size() const {
        return 0;
    }

    bool areNew() const {
        return false;
    }

};

template<typename... Evaluables>
levi::DependenciesSet<0, Evaluables...> levi::make_dependencies_set(levi::ExpressionComponent<Evaluables>&& ...dependencies) {
    return levi::DependenciesSet<0, Evaluables...>(std::forward<levi::ExpressionComponent<Evaluables>>(dependencies)...);
}

template<typename... Evaluables>
levi::DependenciesSet<0, Evaluables...> levi::make_dependencies_set(const levi::ExpressionComponent<Evaluables>& ...dependencies) {
    return levi::DependenciesSet<0, Evaluables...>(dependencies...);
}


#endif // LEVI_DEPENDENCIESSET_H
