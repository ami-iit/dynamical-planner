/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <DynamicalPlanner/Solver.h>
#include <DynamicalPlannerPrivate/Costs.h>
#include <DynamicalPlannerPrivate/Constraints.h>
#include <DynamicalPlannerPrivate/DynamicalConstraints.h>
#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/SharedKinDynComputations.h>

#include <iDynTree/OptimalControlProblem.h>
#include <iDynTree/OCSolvers/MultipleShootingSolver.h>
#include <iDynTree/Integrators/ImplicitTrapezoidal.h>
#include <iDynTree/Optimizers/IpoptInterface.h>

#include <cassert>
#include <iostream>

using namespace DynamicalPlanner;
using namespace DynamicalPlanner::Private;

typedef struct {
    std::shared_ptr<DynamicalConstraints> dynamical;
    std::shared_ptr<CentroidalMomentumConstraint> centroidalMomentum;
    std::shared_ptr<CoMPositionConstraint> comPosition;
    std::vector<std::shared_ptr<ContactVelocityControlConstraints>> leftContactsVelocityControl, rightContactsVelocityControl;
    std::vector<std::shared_ptr<ContactForceControlConstraints>> leftContactsForceControl, rightContactsForceControl;
    std::vector<std::shared_ptr<ContactFrictionConstraint>> leftContactsFriction, rightContactsFriction;
    std::vector<std::shared_ptr<ContactPositionConsistencyConstraint>> leftContactsPosition, rightContactsPosition;
    std::shared_ptr<FeetLateralDistanceConstraint> feetLateralDistance;
    std::shared_ptr<QuaternionNormConstraint> quaternionNorm;
} ConstraintSet;

typedef struct {
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> comPosition;
    std::shared_ptr<FrameOrientationCost> frameOrientation;
    std::vector<std::shared_ptr<ForceMeanCost>> leftForceMeans, rightForceMeans;
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> jointsRegularization;
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> jointsVelocity;
    std::shared_ptr<StaticTorquesCost> staticTorques;
    std::vector<std::shared_ptr<iDynTree::optimalcontrol::L2NormCost>> leftPointsForceDerivative, rightPointsForceDerivative;
    std::vector<std::shared_ptr<iDynTree::optimalcontrol::L2NormCost>> leftPointsAcceleration, rightPointsAcceleration;
} CostsSet;

class Solver::Implementation {
public:
    SettingsStruct settings;

    std::shared_ptr<iDynTree::optimization::Optimizer> optimizer = nullptr;
    std::shared_ptr<iDynTree::optimalcontrol::integrators::Integrator> integrationMethod = nullptr;
    std::shared_ptr<iDynTree::optimalcontrol::OptimalControlProblem> ocProblem;
    std::shared_ptr<iDynTree::optimalcontrol::MultipleShootingSolver> multipleShootingSolver;
    std::shared_ptr<SharedKinDynComputation> sharedKinDyn;

    CostsSet costs;
    ConstraintSet constraints;

    VariablesLabeller stateStructure, controlStructure;

    State initialState, optimalState;
    Control optimalControl;


    bool setVariables(size_t numberOfDofs, size_t numberOfPoints) {

        stateStructure.clear();
        controlStructure.clear();

        setFootVariables("Left", numberOfPoints);
        setFootVariables("Right", numberOfPoints);

        bool ok;
        ok = stateStructure.addLabel("Momentum", 6);
        if (!ok) {
            return false;
        }

        ok = stateStructure.addLabel("CoMPosition", 3);
        if (!ok) {
            return false;
        }

        ok = stateStructure.addLabel("BasePosition", 3);
        if (!ok) {
            return false;
        }

        ok = stateStructure.addLabel("BaseQuaternion", 4);
        if (!ok) {
            return false;
        }

        ok = stateStructure.addLabel("JointsPosition", numberOfDofs);
        if (!ok) {
            return false;
        }

        ok = controlStructure.addLabel("BaseVelocity", 6);
        if (!ok) {
            return false;
        }

        ok = controlStructure.addLabel("JointsVelocity", numberOfDofs);
        if (!ok) {
            return false;
        }

        return true;
    }

    bool setCosts(const SettingsStruct& st, const std::shared_ptr<iDynTree::optimalcontrol::OptimalControlProblem> ocp) {
        bool ok = false;

        if (st.comCostActive) {
            costs.comPosition = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("CoMCost", stateStructure.getIndexRange("CoMPosition"),
                                                                                       stateStructure.size(),iDynTree::IndexRange::InvalidRange(),
                                                                                       controlStructure.size());
            ok = costs.comPosition->setStateWeight(st.comWeights);
            if (!ok) {
                return false;
            }

            ok = costs.comPosition->setStateDesiredTrajectory(st.desiredCoMTrajectory);
            if (!ok) {
                return false;
            }

            ok = ocp->addLagrangeTerm(st.comCostOverallWeight, costs.comPosition);
            if (!ok) {
                return false;
            }
        }

        if (st.frameCostActive) {
            costs.frameOrientation = std::make_shared<FrameOrientationCost>(stateStructure, controlStructure, sharedKinDyn, st.robotModel.getFrameIndex(st.frameForOrientationCost));

            ok = costs.frameOrientation->setDesiredRotationTrajectory(st.desiredRotationTrajectory);
            if (!ok) {
                return false;
            }

            ok = ocp->addLagrangeTerm(st.frameCostOverallWeight, costs.frameOrientation);
            if (!ok) {
                return false;
            }
        }

        if (st.forceMeanCostActive) {
            costs.leftForceMeans.resize(st.leftPointsPosition.size());
            for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {
                costs.leftForceMeans[i] = std::make_shared<ForceMeanCost>(stateStructure, controlStructure, "Left", i);
                ok = ocProblem->addLagrangeTerm(st.forceMeanCostOverallWeight, costs.leftForceMeans[i]);
                if (!ok) {
                    return false;
                }
            }
            costs.rightForceMeans.resize(st.rightPointsPosition.size());
            for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {
                costs.rightForceMeans[i] = std::make_shared<ForceMeanCost>(stateStructure, controlStructure, "Right", i);
                ok = ocProblem->addLagrangeTerm(st.forceMeanCostOverallWeight, costs.rightForceMeans[i]);
                if (!ok) {
                    return false;
                }
            }
        }

        if (st.jointsRegularizationCostActive) {
            costs.jointsRegularization = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("JointsRegularizations", stateStructure.getIndexRange("JointsPosition"),
                                                                                                stateStructure.size(), iDynTree::IndexRange::InvalidRange(),
                                                                                                controlStructure.size());
            ok = costs.jointsRegularization->setStateWeight(st.jointsRegularizationWeights);
            if (!ok) {
                return false;
            }

            ok = costs.jointsRegularization->setStateDesiredTrajectory(st.desiredJointsTrajectory);
            if (!ok) {
                return false;
            }

            ok = ocp->addLagrangeTerm(st.jointsRegularizationCostOverallWeight, costs.jointsRegularization);
            if (!ok) {
                return false;
            }
        }

        if (st.jointsVelocityCostActive) {
            costs.jointsVelocity = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("JointsVelocity", iDynTree::IndexRange::InvalidRange(),
                                                                                          stateStructure.size(), controlStructure.getIndexRange("JointsVelocity"),
                                                                                          controlStructure.size());
            ok = costs.jointsVelocity->setControlWeight(st.jointsVelocityCostWeights);
            if (!ok) {
                return false;
            }

            ok = costs.jointsVelocity->setControlDesiredTrajectory(st.desiredJointsVelocityTrajectory);
            if (!ok) {
                return false;
            }

            ok = ocp->addLagrangeTerm(st.jointsVelocityCostOverallWeight, costs.jointsVelocity);
            if (!ok) {
                return false;
            }
        }

        if (st.staticTorquesCostActive) {
            costs.staticTorques = std::make_shared<StaticTorquesCost>(stateStructure, controlStructure, sharedKinDyn, st.robotModel.getFrameIndex(st.leftFrameName),
                                                                      st.robotModel.getFrameIndex(st.rightFrameName), st.leftPointsPosition, st.rightPointsPosition);

            ok = costs.staticTorques->setWeights(st.staticTorquesCostWeights);
            if (!ok) {
                return false;
            }

            ok = ocp->addLagrangeTerm(st.staticTorquesCostOverallWeight, costs.staticTorques);
            if (!ok) {
                return false;
            }
        }

        if (st.forceDerivativeCostActive) {
            costs.leftPointsForceDerivative.resize(st.leftPointsPosition.size());
            for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {
                costs.leftPointsForceDerivative[i] = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("ForceDerivativeLeft" + std::to_string(i),
                                                                                                            iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                                                            controlStructure.getIndexRange("LeftForceControlPoint" + std::to_string(i)),
                                                                                                            controlStructure.size());
                ok = costs.leftPointsForceDerivative[i]->setControlWeight(st.forceDerivativeWeights);
                if (!ok) {
                    return false;
                }

                ok = costs.leftPointsForceDerivative[i]->setControlDesiredTrajectory(st.desiredForceDerivativeTrajectory);
                if (!ok) {
                    return false;
                }

                ok = ocp->addLagrangeTerm(st.forceDerivativesCostOverallWeight, costs.leftPointsForceDerivative[i]);
                if (!ok) {
                    return false;
                }
            }
            costs.rightPointsForceDerivative.resize(st.rightPointsPosition.size());
            for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {
                costs.rightPointsForceDerivative[i] = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("ForceDerivativeRight" + std::to_string(i),
                                                                                                            iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                                                            controlStructure.getIndexRange("RightForceControlPoint" + std::to_string(i)),
                                                                                                            controlStructure.size());
                ok = costs.rightPointsForceDerivative[i]->setControlWeight(st.forceDerivativeWeights);
                if (!ok) {
                    return false;
                }

                ok = costs.rightPointsForceDerivative[i]->setControlDesiredTrajectory(st.desiredForceDerivativeTrajectory);
                if (!ok) {
                    return false;
                }

                ok = ocp->addLagrangeTerm(st.forceDerivativesCostOverallWeight, costs.rightPointsForceDerivative[i]);
                if (!ok) {
                    return false;
                }
            }
        }

        if (st.pointAccelerationCostActive) {
            costs.leftPointsAcceleration.resize(st.leftPointsPosition.size());
            for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {
                costs.leftPointsAcceleration[i] = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("AccelerationLeft" + std::to_string(i),
                                                                                                            iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                                                            controlStructure.getIndexRange("LeftVelocityControlPoint" + std::to_string(i)),
                                                                                                            controlStructure.size());
                ok = costs.leftPointsAcceleration[i]->setControlWeight(st.pointAccelerationWeights);
                if (!ok) {
                    return false;
                }

                ok = costs.leftPointsAcceleration[i]->setControlDesiredTrajectory(st.desiredPointAccelerationTrajectory);
                if (!ok) {
                    return false;
                }

                ok = ocp->addLagrangeTerm(st.pointAccelerationCostOverallWeight, costs.leftPointsAcceleration[i]);
                if (!ok) {
                    return false;
                }
            }
            costs.rightPointsAcceleration.resize(st.rightPointsPosition.size());
            for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {
                costs.rightPointsAcceleration[i] = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("AccelerationRight" + std::to_string(i),
                                                                                                            iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                                                            controlStructure.getIndexRange("RightVelocityControlPoint" + std::to_string(i)),
                                                                                                            controlStructure.size());
                ok = costs.rightPointsAcceleration[i]->setControlWeight(st.pointAccelerationWeights);
                if (!ok) {
                    return false;
                }

                ok = costs.rightPointsAcceleration[i]->setControlDesiredTrajectory(st.desiredPointAccelerationTrajectory);
                if (!ok) {
                    return false;
                }

                ok = ocp->addLagrangeTerm(st.pointAccelerationCostOverallWeight, costs.rightPointsAcceleration[i]);
                if (!ok) {
                    return false;
                }
            }
        }

        return true;
    }

    bool setConstraints(const SettingsStruct& st, const std::shared_ptr<iDynTree::optimalcontrol::OptimalControlProblem> ocp) {

        HyperbolicSecant forceActivation, velocityActivationXY, velocityActivationZ;
        forceActivation.setScaling(st.forceHyperbolicSecantScaling);
        velocityActivationXY.setScaling(st.velocityHyperbolicSecantScalingXY);
        velocityActivationZ.setScaling(st.velocityHyperbolicSecantScalingZ);

        constraints.leftContactsVelocityControl.resize(st.leftPointsPosition.size());
        constraints.leftContactsForceControl.resize(st.leftPointsPosition.size());
        constraints.leftContactsFriction.resize(st.leftPointsPosition.size());
        constraints.leftContactsPosition.resize(st.leftPointsPosition.size());

        constraints.rightContactsVelocityControl.resize(st.rightPointsPosition.size());
        constraints.rightContactsForceControl.resize(st.rightPointsPosition.size());
        constraints.rightContactsFriction.resize(st.rightPointsPosition.size());
        constraints.rightContactsPosition.resize(st.rightPointsPosition.size());

        iDynTree::FrameIndex leftFrame = st.robotModel.getFrameIndex(st.leftFrameName);
        iDynTree::FrameIndex rightFrame = st.robotModel.getFrameIndex(st.rightFrameName);

        bool ok = false;

        for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {
            constraints.leftContactsVelocityControl[i] = std::make_shared<ContactVelocityControlConstraints>(stateStructure, controlStructure,
                                                                                                             "Left", i, velocityActivationXY,
                                                                                                             velocityActivationZ,
                                                                                                             st.velocityMaximumDerivative,
                                                                                                             st.velocityDissipationRatio);
            ok = ocp->addConstraint(constraints.leftContactsVelocityControl[i]);
            if (!ok) {
                return false;
            }

            constraints.leftContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateStructure, controlStructure, "Left",
                                                                                                       i, forceActivation, st.forceMaximumDerivative,
                                                                                                       st.forceDissipationRatio);
            ok = ocp->addConstraint(constraints.leftContactsForceControl[i]);
            if (!ok) {
                return false;
            }

            constraints.leftContactsFriction[i] = std::make_shared<ContactFrictionConstraint>(stateStructure, controlStructure, "Left", i);

            ok = constraints.leftContactsFriction[i]->setFrictionCoefficient(st.frictionCoefficient);
            if (!ok) {
                return false;
            }

            ok = ocp->addConstraint(constraints.leftContactsFriction[i]);
            if (!ok) {
                return false;
            }

            constraints.leftContactsPosition[i] = std::make_shared<ContactPositionConsistencyConstraint>(stateStructure, controlStructure,
                                                                                                         sharedKinDyn, leftFrame, "Left",
                                                                                                         st.leftPointsPosition[i], i);
            ok = ocp->addConstraint(constraints.leftContactsPosition[i]);
            if (!ok) {
                return false;
            }
        }

        for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {

            constraints.rightContactsVelocityControl[i] = std::make_shared<ContactVelocityControlConstraints>(stateStructure, controlStructure,
                                                                                                              "Right", i, velocityActivationXY,
                                                                                                              velocityActivationZ,
                                                                                                              st.velocityMaximumDerivative,
                                                                                                              st.velocityDissipationRatio);
            ok = ocp->addConstraint(constraints.rightContactsVelocityControl[i]);
            if (!ok) {
                return false;
            }

            constraints.rightContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateStructure, controlStructure, "Right",
                                                                                                        i, forceActivation, st.forceMaximumDerivative,
                                                                                                        st.forceDissipationRatio);
            ok = ocp->addConstraint(constraints.rightContactsForceControl[i]);
            if (!ok) {
                return false;
            }

            constraints.rightContactsFriction[i] = std::make_shared<ContactFrictionConstraint>(stateStructure, controlStructure, "Right", i);

            ok = constraints.rightContactsFriction[i]->setFrictionCoefficient(st.frictionCoefficient);
            if (!ok) {
                return false;
            }

            ok = ocp->addConstraint(constraints.rightContactsFriction[i]);
            if (!ok) {
                return false;
            }

            constraints.rightContactsPosition[i] = std::make_shared<ContactPositionConsistencyConstraint>(stateStructure, controlStructure,
                                                                                                          sharedKinDyn, rightFrame, "Right",
                                                                                                          st.rightPointsPosition[i], i);
            ok = ocp->addConstraint(constraints.rightContactsPosition[i]);
            if (!ok) {
                return false;
            }
        }
        constraints.centroidalMomentum = std::make_shared<CentroidalMomentumConstraint>(stateStructure, controlStructure, sharedKinDyn);
        ok = ocp->addConstraint(constraints.centroidalMomentum);
        if (!ok) {
            return false;
        }

        constraints.comPosition = std::make_shared<CoMPositionConstraint>(stateStructure, controlStructure, sharedKinDyn);
        ok = ocp->addConstraint(constraints.comPosition);
        if (!ok) {
            return false;
        }

        constraints.feetLateralDistance = std::make_shared<FeetLateralDistanceConstraint>(stateStructure, controlStructure, sharedKinDyn,
                                                                                          st.indexOfLateralDirection,
                                                                                          st.robotModel.getFrameIndex(st.referenceFrameNameForFeetDistance),
                                                                                          st.robotModel.getFrameIndex(st.otherFrameNameForFeetDistance));
        ok = constraints.feetLateralDistance->setMinimumDistance(st.minimumFeetDistance);

        ok = ocp->addConstraint(constraints.feetLateralDistance);
        if (!ok) {
            return false;
        }

        constraints.quaternionNorm = std::make_shared<QuaternionNormConstraint>(stateStructure, controlStructure);
        ok = ocp->addConstraint(constraints.quaternionNorm);
        if (!ok) {
            return false;
        }

        return true;
    }

private:

    bool setFootVariables(const std::string& footName, size_t numberOfPoints) {
        bool ok = false;
        for (size_t i = 0; i < numberOfPoints; ++i) {
            ok = stateStructure.addLabel(footName + "ForcePoint" + std::to_string(i), 3);
            if (!ok) {
                return false;
            }
            ok = stateStructure.addLabel(footName + "VelocityPoint" + std::to_string(i), 3);
            if (!ok) {
                return false;
            }
            ok = stateStructure.addLabel(footName + "PositionPoint" + std::to_string(i), 3);
            if (!ok) {
                return false;
            }
            ok = controlStructure.addLabel(footName + "VelocityControlPoint" + std::to_string(i), 3);
            if (!ok) {
                return false;
            }
            ok = controlStructure.addLabel(footName + "ForceControlPoint" + std::to_string(i), 3);
            if (!ok) {
                return false;
            }
        }

        return true;
    }
};



Solver::Solver()
    : m_pimpl(new Implementation)
{
    auto ipoptInterface = std::make_shared<iDynTree::optimization::IpoptInterface>();
    if (ipoptInterface->isAvailable()) {
        m_pimpl->optimizer = ipoptInterface;
    }

}

bool Solver::specifySettings(const Settings &settings)
{
    if (!settings.isValid()) {
        std::cerr << "[ERROR][Solver::specifySettings] The specified settings are not valid." << std::endl;
        return false;
    }

    m_pimpl->settings = settings.getSettings();

    const SettingsStruct& st = m_pimpl->settings;

    bool ok = m_pimpl->setVariables(st.robotModel.getNrOfDOFs(), st.leftPointsPosition.size());

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Error while setting the variables structure." << std::endl;
        return false;
    }

    m_pimpl->sharedKinDyn = std::make_shared<SharedKinDynComputation>();
    ok = m_pimpl->sharedKinDyn->loadRobotModel(st.robotModel);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to configure SharedKinDynComputationsObject." << std::endl;
        return false;
    }

    m_pimpl->constraints.dynamical = std::make_shared<DynamicalConstraints>(m_pimpl->stateStructure,
                                                                            m_pimpl->controlStructure,
                                                                            m_pimpl->sharedKinDyn);

    m_pimpl->ocProblem = std::make_shared<iDynTree::optimalcontrol::OptimalControlProblem>();

    ok =  m_pimpl->ocProblem->setDynamicalSystemConstraint(m_pimpl->constraints.dynamical);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set the dynamical system constraints." << std::endl;
        return false;
    }


    //set costs

    ok = m_pimpl->setCosts(st, m_pimpl->ocProblem);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set the costs." << std::endl;
        return false;
    }


    //set constraints

    ok = m_pimpl->setConstraints(st, m_pimpl->ocProblem);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set constraints." << std::endl;
        return false;
    }



    m_pimpl->multipleShootingSolver = std::make_shared<iDynTree::optimalcontrol::MultipleShootingSolver>(m_pimpl->ocProblem);

    ok = m_pimpl->multipleShootingSolver->setStepSizeBounds(st.minimumDt, st.maximumDt);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set the step size bounds." << std::endl;
        return false;
    }

    ok = m_pimpl->multipleShootingSolver->setControlPeriod(st.controlPeriod);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set the control period." << std::endl;
        return false;
    }

    if (m_pimpl->optimizer) {
        ok = m_pimpl->multipleShootingSolver->setOptimizer(m_pimpl->optimizer);

        if (!ok) {
            std::cerr << "[ERROR][Solver::specifySettings] Failed to set the optimizer." << std::endl;
            return false;
        }
    }

    if (m_pimpl->integrationMethod) {
        if (m_pimpl->integrationMethod->dynamicalSystem().expired()) {
            ok = m_pimpl->multipleShootingSolver->setIntegrator(m_pimpl->integrationMethod);

            if (!ok) {
                std::cerr << "[ERROR][Solver::specifySettings] Failed to set the integration method." << std::endl;
                return false;
            }
        } else {
            std::cerr << "[WARNING][Solver::specifySettings] The specified integrator needs to be recreated. Using the default one." << std::endl;
            auto defaultIntegrator = std::make_shared<iDynTree::optimalcontrol::integrators::ImplicitTrapezoidal>();
            ok = m_pimpl->multipleShootingSolver->setIntegrator(defaultIntegrator);
            if (!ok) {
                std::cerr << "[ERROR][Solver::specifySettings] Failed to set the default integration method." << std::endl;
                return false;
            }
        }
    } else {
        auto defaultIntegrator = std::make_shared<iDynTree::optimalcontrol::integrators::ImplicitTrapezoidal>();
        ok = m_pimpl->multipleShootingSolver->setIntegrator(defaultIntegrator);
        if (!ok) {
            std::cerr << "[ERROR][Solver::specifySettings] Failed to set the default integration method." << std::endl;
            return false;
        }
    }

    return true;
}
