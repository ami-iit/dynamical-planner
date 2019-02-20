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
#include <DynamicalPlannerPrivate/Utilities/levi/QuaternionExpressions.h>
#include <DynamicalPlannerPrivate/Utilities/levi/RelativeQuaternionExpression.h>
#include <cassert>
#include <vector>

namespace DynamicalPlanner {
    namespace Private {
        class RelativeQuaternionEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::RelativeQuaternionEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    ExpressionsServer* m_expressionsServer;
    std::string m_baseName, m_targetName;
    iDynTree::FrameIndex m_baseFrame, m_targetFrame;

    levi::Expression m_relativeJacobian;

public:

    RelativeQuaternionEvaluable(ExpressionsServer* expressionsServer,
                                const std::string& baseFrame,
                                const std::string &targetFrame)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (expressionsServer->jointsPosition(), 4, 1, baseFrame + "_rho_" + targetFrame)
          , m_expressionsServer(expressionsServer)
          , m_baseName(baseFrame)
          , m_targetName(targetFrame)
    {
        const iDynTree::Model& model = expressionsServer->model();
        m_baseFrame = model.getFrameIndex(baseFrame);
        assert(m_baseFrame != iDynTree::FRAME_INVALID_INDEX);
        m_targetFrame = model.getFrameIndex(targetFrame);
        assert(m_targetFrame != iDynTree::FRAME_INVALID_INDEX);

        m_relativeJacobian = expressionsServer->relativeLeftJacobian(baseFrame, targetFrame).block(3, 0, 3, expressionsServer->jointsPosition().rows());
    }

    virtual const LEVI_DEFAULT_MATRIX_TYPE& evaluate() final {

        SharedKinDynComputationsPointer kinDyn = m_expressionsServer->currentKinDyn();

        m_evaluationBuffer = iDynTree::toEigen(kinDyn->getRelativeTransform(m_expressionsServer->currentState(), m_baseFrame, m_targetFrame).getRotation().asQuaternion());

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
RelativeQuaternionEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name()) { //m_expression contains the jointsVariable specified in the constructor

        assert(column == 0);

        levi::unused(column);

        levi::Expression thisQuaternion = m_expressionsServer->relativeQuaternion(m_baseName, m_targetName);

        levi::Expression leftQuaternionMap = 0.5 * G_Expression(thisQuaternion.asVariable()).transpose();

        return leftQuaternionMap * m_relativeJacobian;
    } else {
        return levi::Null(4, variable->dimension());
    }

}


levi::Expression DynamicalPlanner::Private::RelativeQuaternionExpression(ExpressionsServer *expressionsServer, const std::string &baseFrame,
                                                                         const std::string &targetFrame)
{

    return levi::ExpressionComponent<RelativeQuaternionEvaluable>(expressionsServer, baseFrame, targetFrame);
}


