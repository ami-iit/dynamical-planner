/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/DynamicalConstraints.h>
#include <private/SpanUtils.h>
#include <private/QuaternionUtils.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/Utils.h>
#include <cassert>
#include <iostream>
#include <regex>

using namespace DynamicalPlanner::Private;

class DynamicalConstraints::Implementation {
public:
    VariablesLabeller stateVariables;
    VariablesLabeller controlVariables;
    VariablesLabeller dynamics;
    double totalMass;
    size_t leftPoints, rightPoints;
    iDynTree::Vector6 gravityVector;
    iDynTree::Rotation baseRotation;
};


DynamicalConstraints::DynamicalConstraints(VariablesLabeller &stateVariables, VariablesLabeller &controlVariables, double totalMass)
   : iDynTree::optimalcontrol::DynamicalSystem (stateVariables.size(), controlVariables.size())
   , m_pimpl(new Implementation)
{
    m_pimpl->stateVariables = stateVariables;
    m_pimpl->controlVariables = controlVariables;
    m_pimpl->dynamics = m_pimpl->stateVariables;
    m_pimpl->totalMass = totalMass;

    m_pimpl->gravityVector.zero();
    m_pimpl->gravityVector(2) = -9.81;

    m_pimpl->leftPoints = 0;
    m_pimpl->rightPoints = 0;

    for (auto label : stateVariables.listOfLabels()) {
        if (label.find("LeftForcePoint") != std::string::npos) {
            m_pimpl->leftPoints++;
        }

        if (label.find("RightForcePoint") != std::string::npos) {
             m_pimpl->rightPoints++;
        }
    }

    //Check that all corresponding variables are available
    for (size_t i; i < m_pimpl->leftPoints; ++i) {
        if (!stateVariables.getIndexRange("LeftForcePoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftForcePoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("LeftVelocityPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftVelocityPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("LeftPositionPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftPositionPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("LeftForceControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftForceControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("LeftVelocityControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "LeftVelocityControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }
    }

    for (size_t i; i < m_pimpl->rightPoints; ++i) {
        if (!stateVariables.getIndexRange("RightForcePoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightForcePoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("RightVelocityPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightVelocityPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!stateVariables.getIndexRange("RightPositionPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightPositionPoint" + std::to_string(i) << " not available among state variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("RightForceControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightForceControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }

        if (!controlVariables.getIndexRange("RightVelocityControlPoint" + std::to_string(i)).isValid()){
            std::cerr << "[ERROR][DynamicalConstraints::DynamicalConstraints] Variable " << "RightVelocityControlPoint" + std::to_string(i) << " not available among control variables." << std::endl;
            assert(false);
        }
    }

}

DynamicalConstraints::~DynamicalConstraints()
{

}

bool DynamicalConstraints::dynamics(const iDynTree::VectorDynSize &state, double /*time*/, iDynTree::VectorDynSize &stateDynamics)
{
    spanToEigen(m_pimpl->stateVariables.values()) = iDynTree::toEigen(state);
    spanToEigen(m_pimpl->dynamics.values()) = iDynTree::toEigen(state);
    spanToEigen(m_pimpl->controlVariables.values()) = iDynTree::toEigen(controlInput());


    spanToEigen(m_pimpl->dynamics("Momentum")) = m_pimpl->totalMass * iDynTree::toEigen(m_pimpl->gravityVector);

    Eigen::Vector3d distance, appliedForce;

    for (size_t i; i < m_pimpl->leftPoints; ++i) {
        spanToEigen(m_pimpl->dynamics("LeftForcePoint" + std::to_string(i))) = spanToEigen(m_pimpl->controlVariables("LeftForceControlPoint" + std::to_string(i)));
        spanToEigen(m_pimpl->dynamics("LeftVelocityPoint" + std::to_string(i))) = spanToEigen(m_pimpl->controlVariables("LeftVelocityControlPoint" + std::to_string(i)));
        spanToEigen(m_pimpl->dynamics("LeftPositionPoint" + std::to_string(i))) = spanToEigen(m_pimpl->controlVariables("LeftVelocityPoint" + std::to_string(i)));

        spanToEigen(m_pimpl->dynamics("Momentum")).topRows<3>() +=  spanToEigen(m_pimpl->stateVariables("LeftForcePoint" + std::to_string(i)));

        distance = spanToEigen(m_pimpl->stateVariables("LeftPositionPoint" + std::to_string(i))) - spanToEigen(m_pimpl->stateVariables("CoMPosition"));
        appliedForce = spanToEigen(m_pimpl->stateVariables("LeftForcePoint" + std::to_string(i)));
        spanToEigen(m_pimpl->dynamics("Momentum")).bottomRows<3>() += distance.cross(appliedForce);
    }

    for (size_t i; i < m_pimpl->rightPoints; ++i) {
        spanToEigen(m_pimpl->dynamics("RightForcePoint" + std::to_string(i))) = spanToEigen(m_pimpl->controlVariables("RightForceControlPoint" + std::to_string(i)));
        spanToEigen(m_pimpl->dynamics("RightVelocityPoint" + std::to_string(i))) = spanToEigen(m_pimpl->controlVariables("RightVelocityControlPoint" + std::to_string(i)));
        spanToEigen(m_pimpl->dynamics("RightPositionPoint" + std::to_string(i))) = spanToEigen(m_pimpl->controlVariables("RightVelocityPoint" + std::to_string(i)));

        spanToEigen(m_pimpl->dynamics("Momentum")).topRows<3>() +=  spanToEigen(m_pimpl->stateVariables("RightForcePoint" + std::to_string(i)));

        distance = spanToEigen(m_pimpl->dynamics("RightPositionPoint" + std::to_string(i))) - spanToEigen(m_pimpl->stateVariables("CoMPosition"));
        appliedForce = spanToEigen(m_pimpl->stateVariables("RightForcePoint" + std::to_string(i)));
        spanToEigen(m_pimpl->dynamics("Momentum")).bottomRows<3>() += distance.cross(appliedForce);
    }

    spanToEigen(m_pimpl->dynamics("CoMPosition")) = spanToEigen(m_pimpl->stateVariables("Momentum")).topRows<3>()/m_pimpl->totalMass;

    iDynTree::Vector4 normalizedQuaternion;
    iDynTree::toEigen(normalizedQuaternion) = spanToEigen(m_pimpl->stateVariables("BaseQuaternion")).normalized();

    assert(QuaternionBoundsRespected(normalizedQuaternion));

    m_pimpl->baseRotation.fromQuaternion(normalizedQuaternion);

    spanToEigen(m_pimpl->dynamics("BasePosition")) = iDynTree::toEigen(m_pimpl->baseRotation) * spanToEigen(m_pimpl->controlVariables("BaseVelocity")).topRows<3>();

    spanToEigen(m_pimpl->dynamics("BaseQuaternion")) = iDynTree::toEigen(QuaternionLeftTrivializedDerivative(normalizedQuaternion)) * spanToEigen(m_pimpl->controlVariables("BaseVelocity")).bottomRows<3>();

    spanToEigen(m_pimpl->dynamics("JointsPosition")) = spanToEigen(m_pimpl->controlVariables("JointsVelocity"));

    iDynTree::toEigen(stateDynamics) = spanToEigen(m_pimpl->dynamics.values());

    return true;
}

bool DynamicalConstraints::dynamicsStateFirstDerivative(const iDynTree::VectorDynSize &state, double time, iDynTree::MatrixDynSize &dynamicsDerivative)
{

}

bool DynamicalConstraints::dynamicsControlFirstDerivative(const iDynTree::VectorDynSize &state, double time, iDynTree::MatrixDynSize &dynamicsDerivative)
{

}
