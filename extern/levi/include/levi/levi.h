/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef LEVI_H
#define LEVI_H

#if defined (_MSC_VER)
    #ifndef _ENABLE_EXTENDED_ALIGNED_STORAGE
        #ifndef _DISABLE_EXTENDED_ALIGNED_STORAGE
            #define _ENABLE_EXTENDED_ALIGNED_STORAGE
        #endif
    #endif
#endif

#include "Evaluable.h"
#include "Variable.h"
#include "Expression.h"
#include "ExpressionImplementation.h"
#include "DependenciesSet.h"

#endif // LEVI_H
