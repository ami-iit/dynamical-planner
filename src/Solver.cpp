/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <levi/levi.h>
#include <DynamicalPlanner/Solver.h>
#include <DynamicalPlannerPrivate/Costs.h>
#include <DynamicalPlannerPrivate/Constraints.h>
#include <DynamicalPlannerPrivate/Constraints/DynamicalConstraints.h>
#include <DynamicalPlannerPrivate/Utilities/VariablesLabeller.h>
#include <DynamicalPlannerPrivate/Utilities/TimelySharedKinDynComputations.h>
#include <DynamicalPlannerPrivate/Utilities/ExpressionsServer.h>
#include <DynamicalPlannerPrivate/Utilities/QuaternionUtils.h>

#include <iDynTree/OptimalControlProblem.h>
#include <iDynTree/OCSolvers/MultipleShootingSolver.h>
#include <iDynTree/Integrators/ImplicitTrapezoidal.h>
#include <iDynTree/Optimizers/IpoptInterface.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <iDynTree/Core/Span.h>

#include <cassert>
#include <iostream>

using namespace DynamicalPlanner;
using namespace DynamicalPlanner::Private;

typedef struct {
    std::shared_ptr<DynamicalConstraints> dynamical;
    std::shared_ptr<CentroidalMomentumConstraint> centroidalMomentum;
    std::shared_ptr<CoMPositionConstraint> comPosition;
    std::vector<std::shared_ptr<NormalVelocityControlConstraints>> leftNormalVelocityControl, rightNormalVelocityControl;
    std::vector<std::shared_ptr<PlanarVelocityControlConstraints>> leftPlanarVelocityControl, rightPlanarVelocityControl;
    std::vector<std::shared_ptr<ContactForceControlConstraints>> leftContactsForceControl, rightContactsForceControl;
    std::vector<std::shared_ptr<DynamicalComplementarityConstraint>> leftComplementarity, rightComplementarity;
    std::vector<std::shared_ptr<ContactFrictionConstraint>> leftContactsFriction, rightContactsFriction;
    std::vector<std::shared_ptr<ContactPositionConsistencyConstraint>> leftContactsPosition, rightContactsPosition;
    std::shared_ptr<FeetLateralDistanceConstraint> feetLateralDistance;
    std::shared_ptr<QuaternionNormConstraint> quaternionNorm;
} ConstraintSet;

typedef struct {
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> comPosition;
    std::shared_ptr<CoMVelocityCost> comVelocity;
    std::shared_ptr<FrameOrientationCost> frameOrientation;
    std::vector<std::shared_ptr<ForceMeanCost>> leftForceMeans, rightForceMeans;
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> jointsRegularization;
    std::shared_ptr<iDynTree::optimalcontrol::L2NormCost> jointsVelocity;
    std::shared_ptr<StaticTorquesCost> staticTorques;
    std::vector<std::shared_ptr<iDynTree::optimalcontrol::L2NormCost>> leftPointsForceDerivative, rightPointsForceDerivative;
    std::vector<std::shared_ptr<iDynTree::optimalcontrol::L2NormCost>> leftPointsAcceleration, rightPointsAcceleration;
    std::vector<std::shared_ptr<SwingCost>> leftSwings, rightSwings;
    std::vector<std::shared_ptr<PhantomForcesCost>> leftPhantomForces, rightPhantomForces;
    std::shared_ptr<MeanPointPositionCost> meanPositionCost;
} CostsSet;

typedef struct {
    std::vector<iDynTree::IndexRange> positionPoints, forcePoints, velocityControlPoints, forceControlPoints;
} FootRanges;

typedef struct {
    FootRanges left, right;
    iDynTree::IndexRange momentum, comPosition, basePosition;
    iDynTree::IndexRange baseQuaternion, jointsPosition, baseLinearVelocity, baseQuaternionDerivative, jointsVelocity;
} VariablesRanges;

class StateGuesses : public iDynTree::optimalcontrol::TimeVaryingVector {
    std::shared_ptr<TimeVaryingState> m_originalGuesses;
    iDynTree::VectorDynSize m_buffer;
    VariablesRanges m_ranges;

    template<typename Vector>
    void setSegment(iDynTree::IndexRange &range, const Vector &original) {
        iDynTree::toEigen(m_buffer).segment(range.offset, range.size) = iDynTree::toEigen(original);
    }

public:

    StateGuesses(std::shared_ptr<TimeVaryingState> originalGuess, const VariablesRanges &ranges)
        : m_originalGuesses(originalGuess)
        , m_buffer(static_cast<unsigned int>(ranges.left.positionPoints.size() * 2 * 6 + 16 + static_cast<size_t>(ranges.jointsPosition.size)))
        , m_ranges(ranges)
    { }

    ~StateGuesses() override;

    const iDynTree::VectorDynSize &get(double time, bool &isValid) override {

        const State &desiredState = m_originalGuesses->get(time, isValid);

        if (!isValid) {
            std::cerr << "[ERROR][StateGuesses::get] Unable to get a valid state at time " << time << "." << std::endl;
            m_buffer.zero();
            return m_buffer;
        }

        if (!desiredState.checkSize(static_cast<size_t>(m_ranges.jointsPosition.size), m_ranges.left.positionPoints.size())) {
            std::cerr << "[ERROR][StateGuesses::get] Unable to get a state with correct dimensions at time " << time << "." << std::endl;
            isValid = false;
            m_buffer.zero();
            return m_buffer;
        }

        isValid = true;

        for (size_t i = 0; i < m_ranges.left.positionPoints.size(); ++i) {
            setSegment(m_ranges.left.positionPoints[i], desiredState.leftContactPointsState[i].pointPosition);
            setSegment(m_ranges.left.forcePoints[i], desiredState.leftContactPointsState[i].pointForce);
        }

        for (size_t i = 0; i < m_ranges.right.positionPoints.size(); ++i) {
            setSegment(m_ranges.right.positionPoints[i], desiredState.rightContactPointsState[i].pointPosition);
            setSegment(m_ranges.right.forcePoints[i], desiredState.rightContactPointsState[i].pointForce);
        }

        setSegment(m_ranges.momentum, desiredState.momentumInCoM);
        setSegment(m_ranges.comPosition, desiredState.comPosition);
        setSegment(m_ranges.basePosition, desiredState.worldToBaseTransform.getPosition());
        setSegment(m_ranges.baseQuaternion, desiredState.worldToBaseTransform.getRotation().asQuaternion());
        setSegment(m_ranges.jointsPosition, desiredState.jointsConfiguration);

        return m_buffer;
    }
};
StateGuesses::~StateGuesses() { }


class ControlGuesses : public iDynTree::optimalcontrol::TimeVaryingVector {
    std::shared_ptr<TimeVaryingControl> m_originalGuesses;
    iDynTree::VectorDynSize m_buffer;
    VariablesRanges m_ranges;

    template<typename Vector>
    void setSegment(iDynTree::IndexRange &range, const Vector &original) {
        iDynTree::toEigen(m_buffer).segment(range.offset, range.size) = iDynTree::toEigen(original);
    }

public:

    ControlGuesses(std::shared_ptr<TimeVaryingControl> originalGuess, const VariablesRanges &ranges)
        : m_originalGuesses(originalGuess)
        , m_buffer(static_cast<unsigned int>(ranges.left.positionPoints.size() * 12 + 7 + static_cast<size_t>(ranges.jointsPosition.size)))
        , m_ranges(ranges)
    { }

    ~ControlGuesses() override;

    const iDynTree::VectorDynSize &get(double time, bool &isValid) override {

        const Control &desiredControl = m_originalGuesses->get(time, isValid);

        if (!isValid) {
            m_buffer.zero();
            std::cerr << "[ERROR][ControlGuesses::get] Unable to get a valid control at time " << time << "." << std::endl;
            return m_buffer;
        }

        if (!desiredControl.checkSize(static_cast<size_t>(m_ranges.jointsPosition.size), m_ranges.left.positionPoints.size())) {
            std::cerr << "[ERROR][ControlGuesses::get] Unable to get a control with correct dimensions at time " << time << "." << std::endl;
            isValid = false;
            m_buffer.zero();
            return m_buffer;
        }

        isValid = true;

        for (size_t i = 0; i < m_ranges.left.forceControlPoints.size(); ++i) {
            setSegment(m_ranges.left.forceControlPoints[i], desiredControl.leftContactPointsControl[i].pointForceControl);
            setSegment(m_ranges.left.velocityControlPoints[i], desiredControl.leftContactPointsControl[i].pointVelocityControl);
        }

        for (size_t i = 0; i < m_ranges.right.forceControlPoints.size(); ++i) {
            setSegment(m_ranges.right.forceControlPoints[i], desiredControl.rightContactPointsControl[i].pointForceControl);
            setSegment(m_ranges.right.velocityControlPoints[i], desiredControl.rightContactPointsControl[i].pointVelocityControl);
        }

        setSegment(m_ranges.baseLinearVelocity, desiredControl.baseLinearVelocity);
        setSegment(m_ranges.baseQuaternionDerivative, desiredControl.baseQuaternionDerivative);
        setSegment(m_ranges.jointsVelocity, desiredControl.jointsVelocity);

        return m_buffer;
    }
};
ControlGuesses::~ControlGuesses() { }

class VariableBound : public iDynTree::optimalcontrol::TimeVaryingVector {
    iDynTree::VectorDynSize m_firstBounds;
    iDynTree::VectorDynSize m_secondBounds;
    iDynTree::VectorDynSize m_outputBounds;
    double m_switchTime;
    iDynTree::optimalcontrol::TimeRange m_constrainTargetCoMPositionRange;
    iDynTree::IndexRange m_comRange;
    double m_signForTolerance;
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> m_desiredCoMTrajectory;
    std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> m_targetCoMPositionTolerance;
    iDynTree::Vector3 m_toleranceVector;

public:
    VariableBound(iDynTree::VectorDynSize& firstBounds, iDynTree::VectorDynSize& secondBounds, double switchTime)
        : m_firstBounds(firstBounds)
        , m_secondBounds(secondBounds)
        , m_outputBounds(firstBounds)
        , m_switchTime(switchTime)
        , m_desiredCoMTrajectory(nullptr)
    {}

    VariableBound(iDynTree::VectorDynSize& firstBounds, iDynTree::VectorDynSize& secondBounds, double switchTime,
                  const iDynTree::optimalcontrol::TimeRange &constrainTargetCoMPositionRange,
                  const iDynTree::IndexRange& comRange, double signForTolerance,
                  std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingVector> desiredCoMTrajectory,
                  std::shared_ptr<iDynTree::optimalcontrol::TimeVaryingDouble> targetCoMPositionTolerance)
        : m_firstBounds(firstBounds)
        , m_secondBounds(secondBounds)
        , m_outputBounds(firstBounds)
        , m_switchTime(switchTime)
        , m_constrainTargetCoMPositionRange(constrainTargetCoMPositionRange)
        , m_comRange(comRange)
        , m_signForTolerance(signForTolerance)
        , m_desiredCoMTrajectory(desiredCoMTrajectory)
        , m_targetCoMPositionTolerance(targetCoMPositionTolerance)
    {
        assert(m_comRange.isValid());
    }

    virtual ~VariableBound() override;

    virtual const iDynTree::VectorDynSize& get(double time, bool& isValid) override {
        isValid = true;

        if (m_desiredCoMTrajectory && m_targetCoMPositionTolerance && m_constrainTargetCoMPositionRange.isInRange(time)) {
            const iDynTree::VectorDynSize& comReference = m_desiredCoMTrajectory->get(time, isValid);
            if (!isValid) {
                return (time < m_switchTime) ? m_firstBounds : m_secondBounds;
            }

            bool isValidTolerance = false;
            double tolerance = m_targetCoMPositionTolerance->get(time, isValidTolerance);
            if (!isValidTolerance) {
                tolerance = 0;
            }

            m_outputBounds = (time < m_switchTime) ? m_firstBounds : m_secondBounds;

            iDynTree::toEigen(m_toleranceVector).setConstant(m_signForTolerance * tolerance);

            iDynTree::toEigen(m_outputBounds).segment(m_comRange.offset, 3) = iDynTree::toEigen(comReference) + iDynTree::toEigen(m_toleranceVector);

            isValid = true;

            return m_outputBounds;

        } else {
            return (time < m_switchTime) ? m_firstBounds : m_secondBounds;
        }
    }
};
VariableBound::~VariableBound(){}

class Solver::Implementation {
public:
    SettingsStruct settings;

    std::shared_ptr<iDynTree::optimization::Optimizer> optimizer = nullptr;
    std::shared_ptr<iDynTree::optimalcontrol::integrators::Integrator> integrationMethod = nullptr;
    std::shared_ptr<iDynTree::optimalcontrol::OptimalControlProblem> ocProblem;
    std::shared_ptr<iDynTree::optimalcontrol::MultipleShootingSolver> multipleShootingSolver;
    std::shared_ptr<TimelySharedKinDynComputations> timelySharedKinDyn;
    std::shared_ptr<ExpressionsServer> expressionsServer;

    CostsSet costs;
    ConstraintSet constraints;

    VariablesLabeller stateStructure, controlStructure;

    VariablesRanges ranges;
    iDynTree::Position basePositionBuffer;
    iDynTree::Rotation baseRotationBuffer;
    iDynTree::Vector4 baseQuaternion;

    iDynTree::VectorDynSize stateLowerBound, stateUpperBound, controlLowerBound, controlUpperBound, initialStateVector;
    iDynTree::VectorDynSize secondStateLowerBound, secondStateUpperBound, secondControlLowerBound, secondControlUpperBound;
    double plusInfinity, minusInfinity;

    State initialState;
    std::vector<State> optimalStates;
    std::vector<Control> optimalControls;
    std::vector<iDynTree::VectorDynSize> unstructuredOptimalStates;
    std::vector<iDynTree::VectorDynSize> unstructuredOptimalControl;
    std::vector<double> stateTimings;
    std::vector<double> controlTimings;

    std::shared_ptr<StateGuesses> stateGuess;
    std::shared_ptr<ControlGuesses> controlGuess;

    bool prepared;


    bool setVariablesStructure(size_t numberOfDofs, size_t numberOfPoints) {

        stateStructure.clear();
        controlStructure.clear();

        setFootVariables("Left", numberOfPoints, ranges.left);
        setFootVariables("Right", numberOfPoints, ranges.right);

        ranges.momentum = stateStructure.addLabelAndGetIndexRange("Momentum", 6);
        if (!ranges.momentum.isValid()) {
            return false;
        }

        ranges.comPosition = stateStructure.addLabelAndGetIndexRange("CoMPosition", 3);
        if (!ranges.comPosition.isValid()) {
            return false;
        }

        ranges.basePosition = stateStructure.addLabelAndGetIndexRange("BasePosition", 3);
        if (!ranges.basePosition.isValid()) {
            return false;
        }

        ranges.baseQuaternion = stateStructure.addLabelAndGetIndexRange("BaseQuaternion", 4);
        if (!ranges.baseQuaternion.isValid()) {
            return false;
        }

        ranges.jointsPosition = stateStructure.addLabelAndGetIndexRange("JointsPosition", numberOfDofs);
        if (!ranges.jointsPosition.isValid()) {
            return false;
        }

        ranges.baseLinearVelocity = controlStructure.addLabelAndGetIndexRange("BaseLinearVelocity", 3);
        if (!ranges.baseLinearVelocity.isValid()) {
            return false;
        }

        ranges.baseQuaternionDerivative = controlStructure.addLabelAndGetIndexRange("BaseQuaternionDerivative", 4);
        if (!ranges.baseQuaternionDerivative.isValid()) {
            return false;
        }

        ranges.jointsVelocity = controlStructure.addLabelAndGetIndexRange("JointsVelocity", numberOfDofs);
        if (!ranges.jointsVelocity.isValid()) {
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

            ok = ocp->addLagrangeTerm(st.comCostOverallWeight, st.comCostActiveRange, costs.comPosition);
            if (!ok) {
                return false;
            }
        }

        if (st.comVelocityCostActive) {
            costs.comVelocity = std::make_shared<CoMVelocityCost>(stateStructure, controlStructure, stateStructure.getIndexRange("Momentum").offset,
                                                                  timelySharedKinDyn);
            ok = costs.comVelocity->setStateWeight(st.comVelocityWeights);
            if (!ok) {
                return false;
            }

           costs.comVelocity->setLinearVelocityReference(st.desiredCoMVelocityTrajectory);

            ok = ocp->addLagrangeTerm(st.comVelocityCostOverallWeight,
                                      static_cast<std::shared_ptr<iDynTree::optimalcontrol::L2NormCost>>(costs.comVelocity));
            if (!ok) {
                return false;
            }

        }

        if (st.frameCostActive) {
            costs.frameOrientation = std::make_shared<FrameOrientationCost>(stateStructure, controlStructure, timelySharedKinDyn,
                                                                            expressionsServer,
                                                                            st.robotModel.getFrameIndex(st.frameForOrientationCost));

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
            costs.jointsRegularization = std::make_shared<iDynTree::optimalcontrol::L2NormCost>("JointsRegularizations",
                                                                                                stateStructure.getIndexRange("JointsPosition"),
                                                                                                stateStructure.size(),
                                                                                                iDynTree::IndexRange::InvalidRange(),
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
                                                                                          stateStructure.size(),
                                                                                          controlStructure.getIndexRange("JointsVelocity"),
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
            costs.staticTorques = std::make_shared<StaticTorquesCost>(stateStructure, controlStructure, timelySharedKinDyn,
                                                                      st.robotModel.getFrameIndex(st.leftFrameName),
                                                                      st.robotModel.getFrameIndex(st.rightFrameName),
                                                                      st.leftPointsPosition, st.rightPointsPosition);

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
                costs.leftPointsForceDerivative[i] =
                        std::make_shared<iDynTree::optimalcontrol::L2NormCost>("ForceDerivativeLeft" + std::to_string(i),
                                                                               iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                               controlStructure.getIndexRange("LeftForceControlPoint" +
                                                                                                              std::to_string(i)),
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
                costs.rightPointsForceDerivative[i] =
                        std::make_shared<iDynTree::optimalcontrol::L2NormCost>("ForceDerivativeRight" + std::to_string(i),
                                                                               iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                               controlStructure.getIndexRange("RightForceControlPoint" +
                                                                                                              std::to_string(i)),
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
                costs.leftPointsAcceleration[i] =
                        std::make_shared<iDynTree::optimalcontrol::L2NormCost>("AccelerationLeft" + std::to_string(i),
                                                                               iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                               controlStructure.getIndexRange("LeftVelocityControlPoint" +
                                                                                                              std::to_string(i)),
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
                costs.rightPointsAcceleration[i] =
                        std::make_shared<iDynTree::optimalcontrol::L2NormCost>("AccelerationRight" + std::to_string(i),
                                                                               iDynTree::IndexRange::InvalidRange(), stateStructure.size(),
                                                                               controlStructure.getIndexRange("RightVelocityControlPoint" +
                                                                                                              std::to_string(i)),
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

        if (st.swingCostActive) {
            costs.leftSwings.resize(st.leftPointsPosition.size());
            for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {
                costs.leftSwings[i] = std::make_shared<SwingCost>(stateStructure, controlStructure, "Left", i, st.desiredSwingHeight);
                ok = ocProblem->addLagrangeTerm(st.swingCostOverallWeight, costs.leftSwings[i]);
                if (!ok) {
                    return false;
                }
            }
            costs.rightSwings.resize(st.rightPointsPosition.size());
            for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {
                costs.rightSwings[i] = std::make_shared<SwingCost>(stateStructure, controlStructure, "Right", i, st.desiredSwingHeight);
                ok = ocProblem->addLagrangeTerm(st.swingCostOverallWeight, costs.rightSwings[i]);
                if (!ok) {
                    return false;
                }
            }
        }

        if (st.phantomForcesCostActive) {
            HyperbolicSecant forceActivation;
            forceActivation.setScaling(st.normalForceHyperbolicSecantScaling);

            costs.leftPhantomForces.resize(st.leftPointsPosition.size());
            for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {
                costs.leftPhantomForces[i] = std::make_shared<PhantomForcesCost>(stateStructure, controlStructure, "Left", i,
                                                                                 forceActivation);
                ok = ocProblem->addLagrangeTerm(st.phantomForcesCostOverallWeight, costs.leftPhantomForces[i]);
                if (!ok) {
                    return false;
                }
            }
            costs.rightPhantomForces.resize(st.rightPointsPosition.size());
            for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {
                costs.rightPhantomForces[i] = std::make_shared<PhantomForcesCost>(stateStructure, controlStructure, "Right", i,
                                                                                  forceActivation);
                ok = ocProblem->addLagrangeTerm(st.phantomForcesCostOverallWeight, costs.rightPhantomForces[i]);
                if (!ok) {
                    return false;
                }
            }
        }

        if (st.meanPointPositionCostActive) {
            costs.meanPositionCost = std::make_shared<MeanPointPositionCost>(stateStructure, controlStructure);
            costs.meanPositionCost->setDesiredPositionTrajectory(st.desiredMeanPointPosition);
            costs.meanPositionCost->setTimeVaryingWeight(st.meanPointPositionCostTimeVaryingWeight);
            ok = ocProblem->addLagrangeTerm(st.meanPointPositionCostOverallWeight, st.meanPointPositionCostActiveRange, costs.meanPositionCost);
            if (!ok) {
                return false;
            }
        }


        return true;
    }

    bool setConstraints(const SettingsStruct& st, const std::shared_ptr<iDynTree::optimalcontrol::OptimalControlProblem> ocp) {

        HyperbolicSecant forceActivation, velocityActivationZ;
        HyperbolicTangent velocityActivationXY;
        forceActivation.setScaling(st.normalForceHyperbolicSecantScaling);
        velocityActivationXY.setScaling(st.planarVelocityHyperbolicTangentScaling);
        velocityActivationZ.setScaling(st.normalVelocityHyperbolicSecantScaling);

        constraints.leftNormalVelocityControl.resize(st.leftPointsPosition.size());
        constraints.leftPlanarVelocityControl.resize(st.leftPointsPosition.size());
        constraints.leftContactsForceControl.resize(st.leftPointsPosition.size());
        constraints.leftComplementarity.resize(st.leftPointsPosition.size());
        constraints.leftContactsFriction.resize(st.leftPointsPosition.size());
        constraints.leftContactsPosition.resize(st.leftPointsPosition.size());

        constraints.rightNormalVelocityControl.resize(st.rightPointsPosition.size());
        constraints.rightPlanarVelocityControl.resize(st.rightPointsPosition.size());
        constraints.rightContactsForceControl.resize(st.rightPointsPosition.size());
        constraints.rightComplementarity.resize(st.rightPointsPosition.size());
        constraints.rightContactsFriction.resize(st.rightPointsPosition.size());
        constraints.rightContactsPosition.resize(st.rightPointsPosition.size());

        iDynTree::FrameIndex leftFrame = st.robotModel.getFrameIndex(st.leftFrameName);
        iDynTree::FrameIndex rightFrame = st.robotModel.getFrameIndex(st.rightFrameName);

        bool ok = false;

        constraints.centroidalMomentum = std::make_shared<CentroidalMomentumConstraint>(stateStructure, controlStructure,
                                                                                        timelySharedKinDyn, expressionsServer);
        constraints.centroidalMomentum->setEqualityTolerance(st.centroidalMomentumConstraintTolerance);
        ok = ocp->addConstraint(constraints.centroidalMomentum);
        if (!ok) {
            return false;
        }

        constraints.comPosition = std::make_shared<CoMPositionConstraint>(stateStructure, controlStructure,
                                                                          timelySharedKinDyn, expressionsServer);
        constraints.comPosition->setEqualityTolerance(st.comPositionConstraintTolerance);
        ok = ocp->addConstraint(constraints.comPosition);
        if (!ok) {
            return false;
        }

        constraints.feetLateralDistance = std::make_shared<FeetLateralDistanceConstraint>(stateStructure, controlStructure,
                                                                                          timelySharedKinDyn, expressionsServer,
                                                                                          st.indexOfLateralDirection,
                                                                                          st.robotModel.getFrameIndex(
                                                                                              st.referenceFrameNameForFeetDistance),
                                                                                          st.robotModel.getFrameIndex(
                                                                                              st.otherFrameNameForFeetDistance));
        ok = constraints.feetLateralDistance->setMinimumDistance(st.minimumFeetDistance);

        ok = ocp->addConstraint(constraints.feetLateralDistance);
        if (!ok) {
            return false;
        }

        constraints.quaternionNorm = std::make_shared<QuaternionNormConstraint>(stateStructure, controlStructure);
        constraints.quaternionNorm->setEqualityTolerance(st.quaternionModulusConstraintTolerance);
        ok = ocp->addConstraint(constraints.quaternionNorm);
        if (!ok) {
            return false;
        }

        for (size_t i = 0; i < st.leftPointsPosition.size(); ++i) {

//            constraints.leftPlanarVelocityControl[i] = std::make_shared<PlanarVelocityControlConstraints>(stateStructure, controlStructure,
//                                                                                                          "Left", i,velocityActivationXY,
//                                                                                                          st.velocityMaximumDerivative(0),
//                                                                                                          st.velocityMaximumDerivative(1));
//            ok = ocp->addConstraint(constraints.leftPlanarVelocityControl[i]);
//            if (!ok) {
//                return false;
//            }

//            constraints.leftNormalVelocityControl[i] = std::make_shared<NormalVelocityControlConstraints>(stateStructure, controlStructure,
//                                                                                                          "Left", i, velocityActivationZ,
//                                                                                                          st.velocityMaximumDerivative(2));
//            ok = ocp->addConstraint(constraints.leftNormalVelocityControl[i]);
//            if (!ok) {
//                return false;
//            }

            constraints.leftContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateStructure, controlStructure, "Left",
                                                                                                       i, forceActivation,
                                                                                                       st.forceMaximumDerivative(2),
                                                                                                       st.normalForceDissipationRatio,
                                                                                                       st.horizon * st.activeControlPercentage);
            ok = ocp->addConstraint(constraints.leftContactsForceControl[i]);
            if (!ok) {
                return false;
            }

//            constraints.leftComplementarity[i] = std::make_shared<DynamicalComplementarityConstraint>(stateStructure, controlStructure,
//                                                                                                      "Left", i, st.complementarityDissipation);
//            ok = ocp->addConstraint(constraints.leftComplementarity[i]);
//            if (!ok) {
//                return false;
//            }

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
                                                                                                         timelySharedKinDyn, expressionsServer,
                                                                                                         leftFrame, "Left",
                                                                                                         st.leftPointsPosition[i], i);

            constraints.leftContactsPosition[i]->setEqualityTolerance(st.pointPositionConstraintTolerance);

            ok = ocp->addConstraint(constraints.leftContactsPosition[i]);
            if (!ok) {
                return false;
            }
        }

        for (size_t i = 0; i < st.rightPointsPosition.size(); ++i) {

//            constraints.rightPlanarVelocityControl[i] = std::make_shared<PlanarVelocityControlConstraints>(stateStructure, controlStructure,
//                                                                                                          "Right", i,velocityActivationXY,
//                                                                                                          st.velocityMaximumDerivative(0),
//                                                                                                          st.velocityMaximumDerivative(1));
//            ok = ocp->addConstraint(constraints.rightPlanarVelocityControl[i]);
//            if (!ok) {
//                return false;
//            }

            //            constraints.rightNormalVelocityControl[i] = std::make_shared<NormalVelocityControlConstraints>(stateStructure, controlStructure,
            //                                                                                                          "Right", i, velocityActivationZ,
            //                                                                                                          st.velocityMaximumDerivative(2));
            //            ok = ocp->addConstraint(constraints.rightNormalVelocityControl[i]);
            //            if (!ok) {
            //                return false;
            //            }

            constraints.rightContactsForceControl[i] = std::make_shared<ContactForceControlConstraints>(stateStructure, controlStructure, "Right",
                                                                                                        i, forceActivation,
                                                                                                        st.forceMaximumDerivative(2),
                                                                                                        st.normalForceDissipationRatio,
                                                                                                        st.horizon * st.activeControlPercentage);
            ok = ocp->addConstraint(constraints.rightContactsForceControl[i]);
            if (!ok) {
                return false;
            }

//            constraints.rightComplementarity[i] = std::make_shared<DynamicalComplementarityConstraint>(stateStructure, controlStructure,
//                                                                                                      "Right", i, st.complementarityDissipation);
//            ok = ocp->addConstraint(constraints.rightComplementarity[i]);
//            if (!ok) {
//                return false;
//            }

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
                                                                                                          timelySharedKinDyn, expressionsServer,
                                                                                                          rightFrame, "Right",
                                                                                                          st.rightPointsPosition[i], i);

            constraints.rightContactsPosition[i]->setEqualityTolerance(st.pointPositionConstraintTolerance);

            ok = ocp->addConstraint(constraints.rightContactsPosition[i]);
            if (!ok) {
                return false;
            }
        }

        return true;
    }

    bool setBounds(const SettingsStruct& st) {
        stateLowerBound.resize(static_cast<unsigned int>(stateStructure.size()));
        stateUpperBound.resize(static_cast<unsigned int>(stateStructure.size()));
        iDynTree::toEigen(stateLowerBound).setConstant(minusInfinity);
        iDynTree::toEigen(stateUpperBound).setConstant(plusInfinity);

        controlLowerBound.resize(static_cast<unsigned int>(controlStructure.size()));
        controlUpperBound.resize(static_cast<unsigned int>(controlStructure.size()));
        iDynTree::toEigen(controlLowerBound).setConstant(minusInfinity);
        iDynTree::toEigen(controlUpperBound).setConstant(plusInfinity);

        for (size_t i = 0; i < ranges.left.positionPoints.size(); ++i) {
            segment(stateLowerBound, ranges.left.positionPoints[i])(2) = 0.0;
            segment(stateLowerBound, ranges.left.forcePoints[i])(2) = 0.0;
//            iDynTree::toEigen(segment(controlLowerBound, ranges.left.forceControlPoints[i])) = -iDynTree::toEigen(st.forceMaximumDerivative);
//            iDynTree::toEigen(segment(controlUpperBound, ranges.left.forceControlPoints[i])) = iDynTree::toEigen(st.forceMaximumDerivative);
            iDynTree::toEigen(segment(controlLowerBound, ranges.left.velocityControlPoints[i])) = -iDynTree::toEigen(st.velocityMaximumDerivative);
            iDynTree::toEigen(segment(controlUpperBound, ranges.left.velocityControlPoints[i])) = iDynTree::toEigen(st.velocityMaximumDerivative);
        }

        for (size_t i = 0; i < ranges.right.positionPoints.size(); ++i) {
            segment(stateLowerBound, ranges.right.positionPoints[i])(2) = 0.0;
            segment(stateLowerBound, ranges.right.forcePoints[i])(2) = 0.0;
//            iDynTree::toEigen(segment(controlLowerBound, ranges.right.forceControlPoints[i])) = -iDynTree::toEigen(st.forceMaximumDerivative);
//            iDynTree::toEigen(segment(controlUpperBound, ranges.right.forceControlPoints[i])) = iDynTree::toEigen(st.forceMaximumDerivative);
            iDynTree::toEigen(segment(controlLowerBound, ranges.right.velocityControlPoints[i])) = -iDynTree::toEigen(st.velocityMaximumDerivative);
            iDynTree::toEigen(segment(controlUpperBound, ranges.right.velocityControlPoints[i])) = iDynTree::toEigen(st.velocityMaximumDerivative);
        }


        segment(stateLowerBound, ranges.comPosition)(2) = st.minimumCoMHeight;

        iDynTree::toEigen(segment(stateLowerBound, ranges.momentum)).bottomRows<3>().setConstant(-10);
        iDynTree::toEigen(segment(stateUpperBound, ranges.momentum)).bottomRows<3>().setConstant(10);

//        iDynTree::toEigen(segment(stateLowerBound, ranges.baseQuaternion)).setConstant(-1.0);
//        segment(stateLowerBound, ranges.baseQuaternion)(0) = 0.0;
//        iDynTree::toEigen(segment(stateUpperBound, ranges.baseQuaternion)).setConstant(1.0);

        for (size_t j = 0; j < st.jointsLimits.size(); ++j) {
            segment(stateLowerBound, ranges.jointsPosition)(static_cast<long>(j)) = st.jointsLimits[j].first;
            segment(stateUpperBound, ranges.jointsPosition)(static_cast<long>(j)) = st.jointsLimits[j].second;
        }

        for (size_t j = 0; j < st.jointsVelocityLimits.size(); ++j) {
            segment(controlLowerBound, ranges.jointsVelocity)(static_cast<long>(j)) = st.jointsVelocityLimits[j].first;
            segment(controlUpperBound, ranges.jointsVelocity)(static_cast<long>(j)) = st.jointsVelocityLimits[j].second;
        }

        secondStateLowerBound = stateLowerBound;
        secondStateUpperBound = stateUpperBound;
        secondControlLowerBound = controlLowerBound;
        secondControlUpperBound = controlUpperBound;

        for (size_t i = 0; i < ranges.left.positionPoints.size(); ++i) {
            iDynTree::toEigen(segment(secondControlLowerBound, ranges.left.velocityControlPoints[i])).setZero();
            iDynTree::toEigen(segment(secondControlUpperBound, ranges.left.velocityControlPoints[i])).setZero();
        }

        for (size_t i = 0; i < ranges.right.positionPoints.size(); ++i) {
            iDynTree::toEigen(segment(secondControlLowerBound, ranges.right.velocityControlPoints[i])).setZero();
            iDynTree::toEigen(segment(secondControlUpperBound, ranges.right.velocityControlPoints[i])).setZero();
        }

        std::shared_ptr<VariableBound> variableStateLowerBounds;
        std::shared_ptr<VariableBound> variableStateUpperBounds;
        std::shared_ptr<VariableBound> variableControlLowerBounds;
        std::shared_ptr<VariableBound> variableControlUpperBounds;

        if (st.constrainTargetCoMPosition) {
            variableStateLowerBounds = std::make_shared<VariableBound>(stateLowerBound, secondStateLowerBound, st.horizon * st.activeControlPercentage,
                                                                       st.constrainTargetCoMPositionRange, ranges.comPosition, -1.0,
                                                                       st.desiredCoMTrajectory, st.targetCoMPositionTolerance);
            variableStateUpperBounds = std::make_shared<VariableBound>(stateUpperBound, secondStateUpperBound, st.horizon * st.activeControlPercentage,
                                                                       st.constrainTargetCoMPositionRange, ranges.comPosition, +1.0,
                                                                       st.desiredCoMTrajectory, st.targetCoMPositionTolerance);
        } else {
            variableStateLowerBounds = std::make_shared<VariableBound>(stateLowerBound, secondStateLowerBound, st.horizon * st.activeControlPercentage);
            variableStateUpperBounds = std::make_shared<VariableBound>(stateUpperBound, secondStateUpperBound, st.horizon * st.activeControlPercentage);
        }
        variableControlLowerBounds = std::make_shared<VariableBound>(controlLowerBound, secondControlLowerBound, st.horizon * st.activeControlPercentage);
        variableControlUpperBounds = std::make_shared<VariableBound>(controlUpperBound, secondControlUpperBound, st.horizon * st.activeControlPercentage);

        bool ok = ocProblem->setStateBoxConstraints(variableStateLowerBounds, variableStateUpperBounds);
        if (!ok) {
            return false;
        }

        ok = ocProblem->setControlBoxConstraints(variableControlLowerBounds, variableControlUpperBounds);
        if (!ok) {
            return false;
        }

        return true;
    }

    void resizeSolutionVector(size_t numberOfDofs, size_t numberOfPoints) {

        unstructuredOptimalStates.resize(stateTimings.size(), iDynTree::VectorDynSize(static_cast<unsigned int>(stateStructure.size())));
        unstructuredOptimalControl.resize(controlTimings.size(), iDynTree::VectorDynSize(static_cast<unsigned int>(controlStructure.size())));

        optimalStates.resize(stateTimings.size(), State(numberOfDofs, numberOfPoints));
        optimalControls.resize(controlTimings.size(), Control(numberOfDofs, numberOfPoints));
    }


    void fillSolutionVectors() {
        size_t numberOfDofs = settings.robotModel.getNrOfDOFs();
        size_t numberOfPoints = settings.leftPointsPosition.size();

        optimalStates.resize(stateTimings.size(), State(numberOfDofs, numberOfPoints));
        optimalControls.resize(controlTimings.size(), Control(numberOfDofs, numberOfPoints));

        for (size_t i = 0; i < stateTimings.size(); ++i) {
            setStateFromVariables(unstructuredOptimalStates[i], stateTimings[i], optimalStates[i]);
        }

        for (size_t i = 0; i < controlTimings.size(); ++i) {
            setControlFromVariables(unstructuredOptimalControl[i], controlTimings[i], optimalControls[i]);
        }
    }

    void fillInitialState() {
        initialStateVector.resize(static_cast<unsigned int>(ranges.left.positionPoints.size() * 2 * 6 + 16 +
                                                            static_cast<size_t>(ranges.jointsPosition.size)));
        for (size_t i = 0; i < ranges.left.positionPoints.size(); ++i) {
            setSegment(ranges.left.positionPoints[i], initialState.leftContactPointsState[i].pointPosition, initialStateVector);
            setSegment(ranges.left.forcePoints[i], initialState.leftContactPointsState[i].pointForce, initialStateVector);
        }

        for (size_t i = 0; i < ranges.right.positionPoints.size(); ++i) {
            setSegment(ranges.right.positionPoints[i], initialState.rightContactPointsState[i].pointPosition, initialStateVector);
            setSegment(ranges.right.forcePoints[i], initialState.rightContactPointsState[i].pointForce, initialStateVector);
        }

        setSegment(ranges.momentum, initialState.momentumInCoM, initialStateVector);
        setSegment(ranges.comPosition, initialState.comPosition, initialStateVector);
        setSegment(ranges.basePosition, initialState.worldToBaseTransform.getPosition(), initialStateVector);
        setSegment(ranges.baseQuaternion, initialState.worldToBaseTransform.getRotation().asQuaternion(), initialStateVector);
        setSegment(ranges.jointsPosition, initialState.jointsConfiguration, initialStateVector);
    }

private:

    bool setFootVariables(const std::string& footName, size_t numberOfPoints, FootRanges& footRanges) {
        footRanges.forcePoints.resize(numberOfPoints);
        footRanges.positionPoints.resize(numberOfPoints);
        footRanges.velocityControlPoints.resize(numberOfPoints);
        footRanges.forceControlPoints.resize(numberOfPoints);

        for (size_t i = 0; i < numberOfPoints; ++i) {
            footRanges.forcePoints[i] = stateStructure.addLabelAndGetIndexRange(footName + "ForcePoint" + std::to_string(i), 3);
            if (!footRanges.forcePoints[i].isValid()) {
                return false;
            }
            footRanges.positionPoints[i] = stateStructure.addLabelAndGetIndexRange(footName + "PositionPoint" + std::to_string(i), 3);
            if (!footRanges.positionPoints[i].isValid()) {
                return false;
            }
            footRanges.velocityControlPoints[i] = controlStructure.addLabelAndGetIndexRange(footName + "VelocityControlPoint" + std::to_string(i), 3);
            if (!footRanges.velocityControlPoints[i].isValid()) {
                return false;
            }
            footRanges.forceControlPoints[i] = controlStructure.addLabelAndGetIndexRange(footName + "ForceControlPoint" + std::to_string(i), 3);
            if (!footRanges.forceControlPoints[i].isValid()) {
                return false;
            }
        }

        return true;
    }

    iDynTree::Span<const double> segment(const iDynTree::VectorDynSize &fullVector, const iDynTree::IndexRange& indexRange) {
        return iDynTree::make_span(fullVector).subspan(indexRange.offset, indexRange.size);
    }

    iDynTree::Span<double> segment(iDynTree::VectorDynSize &fullVector, const iDynTree::IndexRange& indexRange) {
        return iDynTree::make_span(fullVector).subspan(indexRange.offset, indexRange.size);
    }

    template<typename Vector>
    void setSegment(iDynTree::IndexRange &range, const Vector &original, iDynTree::VectorDynSize& vectorState) {
        iDynTree::toEigen(vectorState).segment(range.offset, range.size) = iDynTree::toEigen(original);
    }

    void setStateFromVariables(const iDynTree::VectorDynSize &unstructured, double time, State &stateToFill) {
        for (size_t i = 0; i < ranges.left.positionPoints.size(); ++i) {
            stateToFill.leftContactPointsState[i].pointPosition = segment(unstructured, ranges.left.positionPoints[i]);
            iDynTree::toEigen(stateToFill.leftContactPointsState[i].pointForce) = iDynTree::toEigen(segment(unstructured, ranges.left.forcePoints[i]));
        }

        for (size_t i = 0; i < ranges.right.positionPoints.size(); ++i) {
            stateToFill.rightContactPointsState[i].pointPosition = segment(unstructured, ranges.right.positionPoints[i]);
            iDynTree::toEigen(stateToFill.rightContactPointsState[i].pointForce) = iDynTree::toEigen(segment(unstructured, ranges.right.forcePoints[i]));
        }

        stateToFill.momentumInCoM = segment(unstructured, ranges.momentum);
        stateToFill.comPosition = segment(unstructured, ranges.comPosition);

        iDynTree::toEigen(basePositionBuffer) = iDynTree::toEigen(segment(unstructured, ranges.basePosition));
        baseQuaternion = segment(unstructured, ranges.baseQuaternion);
        stateToFill.worldToBaseTransform.setPosition(basePositionBuffer);
        baseRotationBuffer.fromQuaternion(NormalizedQuaternion(baseQuaternion));
        stateToFill.worldToBaseTransform.setRotation(baseRotationBuffer);

        stateToFill.jointsConfiguration = segment(unstructured, ranges.jointsPosition);

        stateToFill.time = time;
    }

    void setControlFromVariables(const iDynTree::VectorDynSize &unstructured, double time, Control &controlToFill) {

        for (size_t i = 0; i < ranges.left.forceControlPoints.size(); ++i) {
            controlToFill.leftContactPointsControl[i].pointForceControl = segment(unstructured, ranges.left.forceControlPoints[i]);
            controlToFill.leftContactPointsControl[i].pointVelocityControl = segment(unstructured, ranges.left.velocityControlPoints[i]);
        }

        for (size_t i = 0; i < ranges.right.forceControlPoints.size(); ++i) {
            controlToFill.rightContactPointsControl[i].pointForceControl = segment(unstructured, ranges.right.forceControlPoints[i]);
            controlToFill.rightContactPointsControl[i].pointVelocityControl = segment(unstructured, ranges.right.velocityControlPoints[i]);
        }

        controlToFill.baseLinearVelocity = segment(unstructured, ranges.baseLinearVelocity);
        controlToFill.baseQuaternionDerivative = segment(unstructured, ranges.baseQuaternionDerivative);
        controlToFill.jointsVelocity = segment(unstructured, ranges.jointsVelocity);
        controlToFill.time = time;
    }
};



Solver::Solver()
    : m_pimpl(std::make_unique<Implementation>())
{
    m_pimpl->plusInfinity = 1E19;
    m_pimpl->minusInfinity = -1E19;

    auto ipoptInterface = std::make_shared<iDynTree::optimization::IpoptInterface>();
    if (ipoptInterface->isAvailable()) {
        m_pimpl->optimizer = ipoptInterface;
        m_pimpl->plusInfinity = ipoptInterface->plusInfinity();
        m_pimpl->minusInfinity = ipoptInterface->minusInfinity();
    }

    m_pimpl->prepared = false;

}

Solver::~Solver()
{ }

bool Solver::specifySettings(const Settings &settings)
{
    if (!settings.isValid()) {
        std::cerr << "[ERROR][Solver::specifySettings] The specified settings are not valid." << std::endl;
        return false;
    }

    m_pimpl->settings = settings.getSettings();

    const SettingsStruct& st = m_pimpl->settings;

    size_t numberOfDofs = st.robotModel.getNrOfDOFs();
    size_t numberOfPoints = st.leftPointsPosition.size();

    if (m_pimpl->initialState.checkSize(0,0)) { //the initial state has not been set yet
        m_pimpl->initialState.resize(numberOfDofs, numberOfPoints);
        m_pimpl->initialState.zero();
    } else if (!m_pimpl->initialState.checkSize(numberOfDofs, numberOfPoints)) { //the initial state has already been set but with different dimensions
        std::cerr << "[WARNING][Solver::specifySettings] The initial state was set but with different dimensions. Setting it to zero." << std::endl;
        m_pimpl->initialState.resize(numberOfDofs, numberOfPoints);
        m_pimpl->initialState.zero();
    }

    bool ok = m_pimpl->setVariablesStructure(numberOfDofs, numberOfPoints);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Error while setting the variables structure." << std::endl;
        return false;
    }

    m_pimpl->timelySharedKinDyn = std::make_shared<TimelySharedKinDynComputations>();
    ok = m_pimpl->timelySharedKinDyn->loadRobotModel(st.robotModel);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to configure SharedKinDynComputationsObject." << std::endl;
        return false;
    }

    m_pimpl->expressionsServer = std::make_shared<ExpressionsServer>(m_pimpl->timelySharedKinDyn);

    HyperbolicTangent velocityActivationXY;
    velocityActivationXY.setScaling(st.planarVelocityHyperbolicTangentScaling);

    m_pimpl->constraints.dynamical = std::make_shared<DynamicalConstraints>(m_pimpl->stateStructure,
                                                                            m_pimpl->controlStructure,
                                                                            m_pimpl->timelySharedKinDyn,
                                                                            m_pimpl->expressionsServer,
                                                                            velocityActivationXY);

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

    ok = m_pimpl->ocProblem->setTimeHorizon(m_pimpl->initialState.time, m_pimpl->initialState.time + st.horizon);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set time horizon." << std::endl;
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

    ok = m_pimpl->setBounds(st);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set the variables bounds." << std::endl;
        return false;
    }


    ok = m_pimpl->multipleShootingSolver->getPossibleTimings(m_pimpl->stateTimings, m_pimpl->controlTimings);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to get the possible timings." << std::endl;
        return false;
    }

    ok = m_pimpl->timelySharedKinDyn->setTimings(m_pimpl->stateTimings);

    if (!ok) {
        std::cerr << "[ERROR][Solver::specifySettings] Failed to set the timings to TimelySharedKinDynComputations object." << std::endl;
        return false;
    }

    m_pimpl->resizeSolutionVector(numberOfDofs, numberOfPoints);

    m_pimpl->prepared = true;

    return true;
}

bool Solver::setInitialState(const State &initialState)
{
    if (!(m_pimpl->prepared)) {
        std::cerr << "[ERROR][Solver::setInitialCondition] First you have to specify the settings." << std::endl;
        return false;
    }

    if (!(m_pimpl->initialState.sameSize(initialState))) {
        std::cerr << "[ERROR][Solver::setInitialCondition] The specified initial state dimensions do not match with the specified settings."
                  << std::endl;
        return false;
    }

    m_pimpl->initialState = initialState;

    return true;
}

bool Solver::setGuesses(std::shared_ptr<TimeVaryingState> stateGuesses, std::shared_ptr<TimeVaryingControl> controlGuesses)
{
    if (!stateGuesses) {
        std::cerr << "[ERROR][Solver::setGuesses] The stateGuesses pointer is empty."
                  << std::endl;
        return false;
    }

    if (!controlGuesses) {
        std::cerr << "[ERROR][Solver::setGuesses] The controlGuesses pointer is empty."
                  << std::endl;
        return false;
    }

    m_pimpl->stateGuess = std::make_shared<StateGuesses>(stateGuesses, m_pimpl->ranges);

    m_pimpl->controlGuess = std::make_shared<ControlGuesses>(controlGuesses, m_pimpl->ranges);

    return true;
}

bool Solver::setOptimizer(std::shared_ptr<iDynTree::optimization::Optimizer> optimizer)
{
    if (!optimizer) {
        std::cerr << "[ERROR][Solver::setOptimizer] The optimizer pointer is empty." << std::endl;
        return false;
    }

    if (m_pimpl->prepared){
        if (!(m_pimpl->multipleShootingSolver->setOptimizer(optimizer))) {
            std::cerr << "[ERROR][Solver::setOptimizer] Failed to set the specified optimizer." << std::endl;
            return false;
        }
    }

    m_pimpl->optimizer = optimizer;
    m_pimpl->plusInfinity = optimizer->plusInfinity();
    m_pimpl->minusInfinity = optimizer->minusInfinity();

    return true;
}

bool Solver::setIntegrator(std::shared_ptr<iDynTree::optimalcontrol::integrators::Integrator> integrationMethod)
{
    if (!integrationMethod) {
        std::cerr << "[ERROR][Solver::setIntegrator] The integrationMethod pointer is empty." << std::endl;
        return false;
    }

    if (m_pimpl->prepared){
        if (!(m_pimpl->multipleShootingSolver->setIntegrator(integrationMethod))) {
            std::cerr << "[ERROR][Solver::setIntegrator] Failed to set the specified integration method." << std::endl;
            return false;
        }
    }

    m_pimpl->integrationMethod = integrationMethod;

    return true;
}

bool Solver::solve(std::vector<State> &optimalStates, std::vector<Control> &optimalControls)
{
    if (!(m_pimpl->prepared)) {
        std::cerr << "[ERROR][Solver::solve] First you have to specify the settings." << std::endl;
        return false;
    }

    if (!m_pimpl->initialState.checkSize(m_pimpl->settings.robotModel.getNrOfDOFs(),
                                         m_pimpl->settings.leftPointsPosition.size())) {
        std::cerr <<"[ERROR][Solver::solve] The specified initial dimensions do not match those of the problem." << std::endl;
        return false;
    }

    m_pimpl->fillInitialState();

    bool ok = false;

    ok = m_pimpl->multipleShootingSolver->setInitialState(m_pimpl->initialStateVector);

    if (!ok) {
        std::cerr << "[ERROR][Solver::solve] Failed to set the initial state." << std::endl;
        return false;
    }

    ok = m_pimpl->ocProblem->setTimeHorizon(m_pimpl->initialState.time, m_pimpl->initialState.time + m_pimpl->settings.horizon);

    if (!ok) {
        std::cerr << "[ERROR][Solver::solve] Failed to set the time horizon." << std::endl;
        return false;
    }

    if (m_pimpl->stateGuess && m_pimpl->controlGuess) {
        ok = m_pimpl->multipleShootingSolver->setGuesses(m_pimpl->stateGuess, m_pimpl->controlGuess);
        if (!ok) {
            std::cerr << "[ERROR][Solver::solve] Failed to set guesses." << std::endl;
            return false;
        }
    }

    ok = m_pimpl->multipleShootingSolver->solve();

    if (!ok) {
        std::cerr << "[ERROR][Solver::solve] Failed to solve the optimization problem." << std::endl;
        return false;
    }

    ok = m_pimpl->multipleShootingSolver->getTimings(m_pimpl->stateTimings, m_pimpl->controlTimings);

    if (!ok) {
        std::cerr << "[ERROR][Solver::solve] Failed to retrieve the timings." << std::endl;
        return false;
    }

    ok = m_pimpl->multipleShootingSolver->getSolution(m_pimpl->unstructuredOptimalStates,
                                                      m_pimpl->unstructuredOptimalControl);

    if (!ok) {
        std::cerr << "[ERROR][Solver::solve] Failed to retrieve the optimal solution." << std::endl;
        return false;
    }

    m_pimpl->fillSolutionVectors();

    optimalStates = m_pimpl->optimalStates;

    optimalControls = m_pimpl->optimalControls;

    m_pimpl->stateGuess = nullptr;
    m_pimpl->controlGuess = nullptr;


    return true;
}

const std::vector<State> &Solver::optimalStates() const
{
    return m_pimpl->optimalStates;
}

const std::vector<Control> &Solver::optimalControls() const
{
    return m_pimpl->optimalControls;
}
