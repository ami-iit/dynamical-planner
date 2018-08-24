/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */
#ifndef DPLANNER_CHECKEQUALVECTOR_H
#define DPLANNER_CHECKEQUALVECTOR_H

#include <iDynTree/Core/Utils.h>

namespace DynamicalPlanner {
    namespace Private {

        template<typename VectorType>
        bool VectorsAreEqual(const VectorType& vec1, const VectorType& vec2, double tol = iDynTree::DEFAULT_TOL){
            if (vec1.size() != vec2.size()) {
                return false;
            }

            for (unsigned int i = 0; i < vec1.size(); ++i) {
                if (!iDynTree::checkDoublesAreEqual(vec1(i), vec2(i), tol))
                    return false;
            }

            return true;
        }

    }
}

#endif // DPLANNER_CHECKEQUALVECTOR_H
