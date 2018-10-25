/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlanner/RectangularFoot.h>
#include <iDynTree/Core/EigenHelpers.h>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace DynamicalPlanner;

class RectangularFoot::Implementation {
public:
    std::vector<iDynTree::Position> pointsPositions;
    iDynTree::Position centerPosition;
    iDynTree::MatrixFixSize<3, 8> regressionMatrix;
    iDynTree::MatrixFixSize<3, 4> reducedRegressionMatrix;
    iDynTree::iDynTreeEigenMatrix regressionMatrixBuffer, reducedRegressionMatrixBuffer;
    Eigen::VectorXd regressionVector;
    iDynTree::VectorFixSize<8> solutionVector;
    iDynTree::VectorFixSize<4> reducedSolutionVector;
    std::vector<bool> pointActive;
    double halfL, halfD;
    bool footSet;
};

RectangularFoot::RectangularFoot()
    : m_pimpl(new Implementation)
{
    m_pimpl->pointsPositions.resize(4);
    m_pimpl->footSet = false;

    m_pimpl->regressionMatrixBuffer.resize(3, 8);
    m_pimpl->reducedRegressionMatrixBuffer.resize(3, 4);

    iDynTree::toEigen(m_pimpl->reducedRegressionMatrix).block<2,2>(0,0).setIdentity();
    iDynTree::toEigen(m_pimpl->reducedRegressionMatrix).block<2,2>(0,2).setIdentity();

    m_pimpl->regressionVector.resize(3);
    m_pimpl->regressionVector.setZero();
    m_pimpl->pointActive.resize(4, true);
}

RectangularFoot::~RectangularFoot()
{ }

bool RectangularFoot::setFoot(double xLength, double yLength, const iDynTree::Position &topLeftPointPosition)
{
    if (xLength <= 0) {
        std::cerr << "[ERROR][RectangularFoot::setFoot] xLength is supposed to be positive." << std::endl;
        return false;
    }

    if (yLength <= 0) {
        std::cerr << "[ERROR][RectangularFoot::setFoot] yLength is supposed to be positive." << std::endl;
        return false;
    }

    m_pimpl->halfL = xLength/2.0;
    m_pimpl->halfD = yLength/2.0;

    iDynTree::Position centerInTopLeftCoordinates(-(m_pimpl->halfL), -(m_pimpl->halfD), 0.0);
    m_pimpl->centerPosition = topLeftPointPosition + centerInTopLeftCoordinates;

    m_pimpl->pointsPositions[0](0) = m_pimpl->halfL;
    m_pimpl->pointsPositions[0](1) = m_pimpl->halfD;
    m_pimpl->pointsPositions[0](2) = 0.0;

    m_pimpl->pointsPositions[1](0) = m_pimpl->halfL;
    m_pimpl->pointsPositions[1](1) = -(m_pimpl->halfD);
    m_pimpl->pointsPositions[1](2) = 0.0;

    m_pimpl->pointsPositions[2](0) = -(m_pimpl->halfL);
    m_pimpl->pointsPositions[2](1) = m_pimpl->halfD;
    m_pimpl->pointsPositions[2](2) = 0.0;

    m_pimpl->pointsPositions[3](0) = -(m_pimpl->halfL);
    m_pimpl->pointsPositions[3](1) = -(m_pimpl->halfD);
    m_pimpl->pointsPositions[3](2) = 0.0;

    for (unsigned int i = 0; i < 4; ++i) {
        iDynTree::toEigen(m_pimpl->regressionMatrix).block<2,2>(0, 2 * i).setIdentity();
        m_pimpl->regressionMatrix(2, 2 * i) = - m_pimpl->pointsPositions[i](1);
        m_pimpl->regressionMatrix(2, 2 * i + 1) = m_pimpl->pointsPositions[i](0);
    }

    m_pimpl->footSet = true;

    return true;
}

bool RectangularFoot::getPoints(const iDynTree::Transform &frameTransform, std::vector<iDynTree::Position> &pointsPosition) const
{
    if (!(m_pimpl->footSet)) {
        std::cerr << "[ERROR][RectangularFoot::getPoints] First you have to call the setFoot method." << std::endl;
        return false;
    }

    pointsPosition.resize(4);

    for (size_t i = 0; i < 4; i++) {
        pointsPosition[i] = frameTransform * (m_pimpl->centerPosition + m_pimpl->pointsPositions[i]);
    }

    return true;
}

bool RectangularFoot::getForces(const iDynTree::Wrench &footWrench, std::vector<iDynTree::Force> &pointsForces)
{
    // See https://github.com/loc2/component_wholebody-teleoperation/issues/107

    if (!(m_pimpl->footSet)) {
        std::cerr << "[ERROR][RectangularFoot::getForces] First you have to call the setFoot method." << std::endl;
        return false;
    }

    pointsForces.resize(4);

    iDynTree::Transform centerToFrameTransform(iDynTree::Rotation::Identity(), -m_pimpl->centerPosition);

    iDynTree::Wrench wrenchInCenter = centerToFrameTransform * footWrench;

    double tx = wrenchInCenter.getAngularVec3()(0);
    double ty = wrenchInCenter.getAngularVec3()(1);
    double tz = wrenchInCenter.getAngularVec3()(2);
    double fx = wrenchInCenter.getLinearVec3()(0);
    double fy = wrenchInCenter.getLinearVec3()(1);
    double fz = wrenchInCenter.getLinearVec3()(2);

    double copX = - ty/fz;
    double copY = tx/fz;

    double l = 2.0 * m_pimpl->halfL;
    double d = 2.0 * m_pimpl->halfD;

    if ((std::fabs(copX) > m_pimpl->halfL) || (std::fabs(copY) > m_pimpl->halfD)) {
        std::cerr << "[WARNING][RectangularFoot::getForces] The COP is outside the foot boundaries. Some of the normal forces may result negative." << std::endl;
    }

    if (fz < 1e-7) {
        for (unsigned int i = 0; i < 4; ++i) {
            pointsForces[i].zero();
        }
        return true;
    }

    iDynTree::Vector4 alphas;

    alphas(3) = (std::max(0.0, -copX/l - copY/d) + std::min(0.5 - copY/d, 0.5 - copX/l))/2;

    pointsForces[3](2) = alphas(3) * fz;
    pointsForces[0](2) = pointsForces[3](2) - ty/l + tx/d;
    alphas(0) = pointsForces[0](2) / fz;
    pointsForces[1](2) = -(pointsForces[3](2)) + fz/2.0 - tx/d;
    alphas(1) = pointsForces[1](2) / fz;
    alphas(2) = 1.0 - alphas(0) - alphas(1) - alphas(3);
    pointsForces[2](2) = fz - pointsForces[3](2) - pointsForces[0](2) - pointsForces[1](2);

    //Having obtained the normal forces we can obtain the remaining components using pseudoinverse

    m_pimpl->regressionVector(0) = fx;
    m_pimpl->regressionVector(1) = fy;
    m_pimpl->regressionVector(2) = tz;

    size_t activePoints = 0;
    for (unsigned int i = 0; i < 4; ++i) {
        if (std::fabs(alphas(i)) > 1e-7) {
            m_pimpl->pointActive[i] = true;
            activePoints++;
        } else {
            m_pimpl->pointActive[i] = false;
        }
    }

    if (activePoints == 1) {

        for (unsigned int i = 0; i < 4; ++i) {
            if (m_pimpl->pointActive[i]) {
                pointsForces[i](0) = fx;
                pointsForces[i](1) = fy;
            } else {
                pointsForces[i].zero();
            }
        }

        return true;
    }

    if (activePoints == 2) {
        unsigned int el = 0;
        for (unsigned int i = 0; i < 4; ++i) {
            if (m_pimpl->pointActive[i]) {
                m_pimpl->reducedRegressionMatrix(2, 2 * el) = - m_pimpl->pointsPositions[i](1);
                m_pimpl->reducedRegressionMatrix(2, 2 * el + 1) = m_pimpl->pointsPositions[i](0);
                el++;
            }
        }
        m_pimpl->reducedRegressionMatrixBuffer = iDynTree::toEigen(m_pimpl->reducedRegressionMatrix);
        iDynTree::toEigen(m_pimpl->reducedSolutionVector) = m_pimpl->reducedRegressionMatrixBuffer.colPivHouseholderQr().solve(m_pimpl->regressionVector);

        el = 0;

        for (unsigned int i = 0; i < 4; ++i) {
            if (m_pimpl->pointActive[i]) {
                pointsForces[i](0) = m_pimpl->reducedSolutionVector(2*el);
                pointsForces[i](1) = m_pimpl->reducedSolutionVector(2*el + 1);
                el++;
            } else {
                pointsForces[i].zero();
            }
        }
        return true;

    }

    m_pimpl->regressionMatrixBuffer = iDynTree::toEigen(m_pimpl->regressionMatrix);
    m_pimpl->solutionVector.zero();
    iDynTree::toEigen(m_pimpl->solutionVector) = m_pimpl->regressionMatrixBuffer.colPivHouseholderQr().solve(m_pimpl->regressionVector);

    for (unsigned int i = 0; i < 4; ++i) {
        pointsForces[i](0) = m_pimpl->solutionVector(2*i);
        pointsForces[i](1) = m_pimpl->solutionVector(2*i + 1);
    }

    return true;
}
