/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <iDynTree/Model/Model.h>
#include <iDynTree/Model/Traversal.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativePositionExpression.h>
#include <cassert>
#include <vector>

namespace DynamicalPlanner {
    namespace Private {
        class RelativePositionEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::RelativePositionEvaluable : public levi::DefaultEvaluable {

    ExpressionsServer* m_expressionsServer;
    levi::ScalarVariable m_timeVariable;
    std::string m_baseName, m_targetName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;

    levi::Expression m_derivative;
    levi::Variable m_jointsVariable;

public:

    RelativePositionEvaluable(ExpressionsServer* expressionsServer,
                              const std::string& baseFrame,
                              const std::string &targetFrame)
        : levi::DefaultEvaluable(3, 1, baseFrame + "_p_" + targetFrame)
          , m_expressionsServer(expressionsServer)
          , m_baseName(baseFrame)
          , m_targetName(targetFrame)
          , m_jointsVariable(expressionsServer->jointsPosition())
    {
        const iDynTree::Model& model = expressionsServer->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);

        m_derivative = expressionsServer->relativeRotation(baseFrame, targetFrame) *
            expressionsServer->relativeLeftJacobian(baseFrame, targetFrame).block(0, 0, 3, expressionsServer->jointsPosition().rows());

        addDependencies(m_jointsVariable);
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        SharedKinDynComputationsPointer kinDyn = m_expressionsServer->currentKinDyn();

        m_evaluationBuffer = iDynTree::toEigen(kinDyn->getRelativeTransform(m_expressionsServer->currentState(), m_baseFrame, m_targetFrame).getPosition());

        return m_evaluationBuffer;
    }

    virtual void clearDerivativesCache() final {
        this->m_derivativeBuffer.clear();
        m_derivative.clearDerivativesCache();
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
RelativePositionEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_jointsVariable.name()) { //m_expression contains the jointsVariable specified in the constructor

        assert(column == 0);

        levi::unused(column);

        return m_derivative;
    } else {
        return levi::Null(3, variable->dimension());
    }

}


levi::Expression DynamicalPlanner::Private::RelativePositionExpression(ExpressionsServer *expressionsServer, const std::string &baseFrame,
                                                                       const std::string &targetFrame)
{

    return levi::ExpressionComponent<RelativePositionEvaluable>(expressionsServer, baseFrame, targetFrame);
}


