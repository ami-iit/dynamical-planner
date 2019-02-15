/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/
#ifndef DPLANNER_TRANSFORMEXPRESSION_H
#define DPLANNER_TRANSFORMEXPRESSION_H

#include <levi/ForwardDeclarations.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <memory>
#include <string>

namespace DynamicalPlanner {
    namespace Private {

        class TransformExpression {

            class Implementation;
            std::unique_ptr<Implementation> m_pimpl;

        public:

            TransformExpression();

            TransformExpression(const levi::Expression& position, const levi::Expression& rotation);

            TransformExpression(const TransformExpression& other);

            ~TransformExpression();

            levi::Expression position() const;

            levi::Expression rotation() const;

            levi::Expression operator*(const levi::Expression& position) const;

            void operator=(const TransformExpression& other);

            static TransformExpression RelativeTransform(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                                                         RobotState *robotState,
                                                         const std::string& baseFrame,
                                                         const std::string &targetFrame,
                                                         levi::Variable jointsVariable,
                                                         levi::ScalarVariable timeVariable);
        };
    }
}

DynamicalPlanner::Private::TransformExpression operator*(const DynamicalPlanner::Private::TransformExpression& lhs, const DynamicalPlanner::Private::TransformExpression& rhs);


#endif // DPLANNER_TRANSFORMEXPRESSION_H
