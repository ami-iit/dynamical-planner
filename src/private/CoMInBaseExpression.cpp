/*
* Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
* Authors: Stefano Dafarra
* CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
*
*/

#include <levi/levi.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <DynamicalPlannerPrivate/Utilities/levi/AdjointTransformExpression.h>
#include <DynamicalPlannerPrivate/Utilities/levi/CoMInBaseExpression.h>

namespace DynamicalPlanner {
    namespace Private {
        class CoMInBaseJacobianEvaluable;
        class CoMInBasePositionEvaluable;
    }
}

using namespace DynamicalPlanner::Private;

class DynamicalPlanner::Private::CoMInBaseJacobianEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    std::shared_ptr<TimelySharedKinDynComputations> m_timelySharedKinDyn;
    RobotState* m_robotState;
    levi::ScalarVariable m_timeVariable;
    iDynTree::MatrixDynSize m_jacobian;


    std::vector<levi::Expression> m_columns;

public:

    CoMInBaseJacobianEvaluable(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                         RobotState *robotState,
                         levi::Variable jointsVariable,
                         levi::ScalarVariable timeVariable)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (jointsVariable, 3, jointsVariable.rows(), "B_J_[B, CoM]")
          , m_timelySharedKinDyn(sharedKinDyn)
          , m_robotState(robotState)
          , m_timeVariable(timeVariable)
    {
        const iDynTree::Model& model = sharedKinDyn->model();

        iDynTree::Traversal traversal;
        iDynTree::LinkIndex baseLink = model.getLinkIndex(sharedKinDyn->getFloatingBase());
        assert(baseLink != iDynTree::LINK_INVALID_INDEX);
        bool ok = model.computeFullTreeTraversal(traversal, baseLink);
        assert(ok);

        double totalMass = 0.0;

        for(size_t l=0; l < model.getNrOfLinks(); l++)
        {
            totalMass += model.getLink(static_cast<iDynTree::LinkIndex>(l))->getInertia().getMass();
        }

        m_columns.resize(static_cast<size_t>(jointsVariable.rows()), levi::Null(3,1));

        iDynTree::IJointConstPtr jointPtr;
        iDynTree::LinkIndex linkIndex;
        iDynTree::LinkConstPtr linkPtr;
        size_t jointIndex;
        std::string linkName;

        iDynTree::LinkIndex visitedLink, childLink, parentLink;

        Eigen::Matrix<double, 6,1> motionSubSpaceVector;

        for (size_t l = 0; l < model.getNrOfLinks(); ++l) {

            linkIndex = static_cast<iDynTree::LinkIndex>(l);
            assert(model.isValidLinkIndex(linkIndex));
            linkName = model.getLinkName(linkIndex);
            linkPtr = model.getLink(linkIndex);
            double relativeMass = linkPtr->getInertia().getMass() / totalMass;

            visitedLink = linkIndex;

            while (visitedLink != baseLink) {
                jointPtr = traversal.getParentJointFromLinkIndex(visitedLink);
                jointIndex = static_cast<size_t>(jointPtr->getIndex());

                childLink =  traversal.getChildLinkIndexFromJointIndex(model, jointPtr->getIndex());
                parentLink = traversal.getParentLinkIndexFromJointIndex(model, jointPtr->getIndex());

                motionSubSpaceVector = iDynTree::toEigen(jointPtr->getMotionSubspaceVector(0,
                                                                                           childLink,
                                                                                           parentLink)) * relativeMass;

                m_columns[jointIndex] = m_columns[jointIndex] + AdjointTransformExpression(sharedKinDyn, robotState, linkName, model.getLinkName(childLink), jointsVariable, timeVariable).block(0,0,3,6) *
                        levi::Constant(motionSubSpaceVector, "s_" + std::to_string(jointIndex) + " m_" + std::to_string(jointIndex) + "/M");

                visitedLink = traversal.getParentLinkFromLinkIndex(visitedLink)->getIndex();
            }
        }

    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_timelySharedKinDyn->get(m_timeVariable.evaluate());


        iDynTree::Vector4 normalizedQuaternion;
        iDynTree::toEigen(normalizedQuaternion) = iDynTree::toEigen(m_robotState->base_quaternion).normalized();
        iDynTree::Rotation baseRotation;
        baseRotation.fromQuaternion(normalizedQuaternion);

        bool ok = kinDyn->getCenterOfMassJacobian(*m_robotState, m_jacobian);
        assert(ok);

        m_evaluationBuffer = iDynTree::toEigen(baseRotation.inverse()) * iDynTree::toEigen(m_jacobian).rightCols(m_robotState->s.size());

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
CoMInBaseJacobianEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name()) { //m_expression contains the jointsVariable specified in the constructor

        return (m_columns[static_cast<size_t>(column)]).getColumnDerivative(0, variable);

    } else {
        return levi::Null(6, variable->dimension());
    }

}

class DynamicalPlanner::Private::CoMInBasePositionEvaluable :
    public levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable> { //Using UnaryOperator as base class, since the transform only depends on joints values. This allows to reuse some buffer mechanism for derivatives and the definition of isNew()

    std::shared_ptr<TimelySharedKinDynComputations> m_timelySharedKinDyn;
    RobotState* m_robotState;
    levi::ScalarVariable m_timeVariable;

    levi::Expression m_derivative;

public:

    CoMInBasePositionEvaluable(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn,
                              RobotState *robotState,
                              levi::Variable jointsVariable,
                              levi::ScalarVariable timeVariable)
        : levi::UnaryOperator<LEVI_DEFAULT_MATRIX_TYPE, levi::DefaultVariableEvaluable>
          (jointsVariable, 3, 1, "b_p_CoM")
          , m_timelySharedKinDyn(sharedKinDyn)
          , m_robotState(robotState)
          , m_timeVariable(timeVariable)
    {

        m_derivative = levi::ExpressionComponent<CoMInBaseJacobianEvaluable>(sharedKinDyn, robotState, jointsVariable, timeVariable);
    }

    const LEVI_DEFAULT_MATRIX_TYPE& evaluate() {

        SharedKinDynComputationsPointer kinDyn = m_timelySharedKinDyn->get(m_timeVariable.evaluate());


        iDynTree::Vector4 normalizedQuaternion;
        iDynTree::toEigen(normalizedQuaternion) = iDynTree::toEigen(m_robotState->base_quaternion).normalized();
        iDynTree::Rotation baseRotation;
        baseRotation.fromQuaternion(normalizedQuaternion);

        iDynTree::Transform baseTransform;
        baseTransform.setRotation(baseRotation);
        iDynTree::Position basePosition;
        iDynTree::toEigen(basePosition) = iDynTree::toEigen(m_robotState->base_position);
        baseTransform.setPosition(basePosition);

        iDynTree::Position comPosition = kinDyn->getCenterOfMassPosition(*m_robotState);

        m_evaluationBuffer = iDynTree::toEigen(baseTransform.inverse() * comPosition);

        return m_evaluationBuffer;
    }

    virtual levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
    getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable) final;
};

levi::ExpressionComponent<typename levi::DefaultEvaluable::derivative_evaluable>
CoMInBasePositionEvaluable::getNewColumnDerivative(Eigen::Index column, std::shared_ptr<levi::VariableBase> variable)
{
    if (variable->variableName() == m_expression.name()) { //m_expression contains the jointsVariable specified in the constructor

        assert(column == 0);

        levi::unused(column);

        return m_derivative;
    } else {
        return levi::Null(3, variable->dimension());
    }

}

levi::Expression DynamicalPlanner::Private::CoMInBaseExpression(std::shared_ptr<TimelySharedKinDynComputations> sharedKinDyn, RobotState *robotState,
                                                                levi::Variable jointsVariable, levi::ScalarVariable timeVariable)
{
    return levi::ExpressionComponent<CoMInBasePositionEvaluable>(sharedKinDyn, robotState, jointsVariable, timeVariable);
}

