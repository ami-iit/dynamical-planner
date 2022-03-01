/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_VARIABLEBASE_H
#define LEVI_VARIABLEBASE_H

#include <levi/HelpersForwardDeclarations.h>
#include <levi/ForwardDeclarations.h>

/**
 * @brief The VariableBase class
 *
 * Base class for variables.
 */
class levi::VariableBase {

protected:
    /**
     * @brief Name of the variable
     */
    std::string m_name;

    /**
     * @brief Dimension of the variable
     */
    Eigen::Index m_dimension;

    /**
     * @brief Constructor
     * @param dimension of the variable
     * @param name of the variable
     *
     * @Note this constructor can be called only by derived classes
     */
    VariableBase(Eigen::Index dimension, const std::string& name)
        : m_name(name)
        , m_dimension(dimension)
    { }

public:

    /**
     * @brief Return the dimension of the variable
     * @return The dimension of the variable.
     */
    Eigen::Index dimension() const{
        return m_dimension;
    }

    /**
     * @brief Return the name of the variable
     * @return The name of the variable.
     */
    std::string variableName() const {
        return m_name;
    }

};

#endif // LEVI_VARIABLEBASE_H
