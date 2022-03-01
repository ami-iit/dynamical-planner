/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <levi/levi.h>
#include <iostream>


int main() {
    using namespace levi;

    Identity identity(3,3);

    Expression antiIdentity = Expression::ComposeByRows({identity.row(2), identity.row(1), identity.row(0)}, "Anti-identity");

    std::cerr << antiIdentity.evaluate() << std::endl;

    Variable x(2, "x");
    auto derivative = antiIdentity.getColumnDerivative(1, x);

    Expression antiIdentity2 = Expression::ComposeByCols({identity.col(2), identity.col(1), identity.col(0)}, "Anti-identity2");

    std::cerr << antiIdentity2.evaluate() << std::endl;

    auto derivative2 = antiIdentity2.getColumnDerivative(1, x);

    Variable k(3, "k");

    Expression kSquared = k.transpose() * k;

    Eigen::Vector3d k_values;
    k_values = Eigen::Vector3d::Random();
    k = k_values;

    Eigen::MatrixXd squaredNorm(1,1);
    squaredNorm(0,0) = k_values.squaredNorm();

    assert((kSquared.evaluate() - squaredNorm).cwiseAbs().maxCoeff() < 1e-10);

    Expression kSquaredDerivative = kSquared.getColumnDerivative(0, k);

    assert((kSquaredDerivative.evaluate() - 2*k_values.transpose()).cwiseAbs().maxCoeff() < 1e-10);

    Expression twiceX;

    twiceX = 2.0 * x;

    Variable x2 = twiceX.asVariable();

    Expression test = x2.transpose() * x2;

    Expression testDerivative = test.getColumnDerivative(0, x2);

    assert(testDerivative.evaluate() == (2.0 * twiceX).evaluate().transpose());

    Expression testDoubleDerivative = testDerivative.getColumnDerivative(0, x);

    assert(testDoubleDerivative.evaluate()(0,0) == 4.0 && testDoubleDerivative.evaluate()(0,1) == 0.0);

    Variable y(3, "y"), w(3, "w");

    y = Eigen::VectorXd::Random(3,1);
    w = Eigen::VectorXd::Random(3,1);

    levi::Expression horzcat = levi::Expression::Horzcat(y.skew(), w.skew(), "horzcat");

    Eigen::MatrixXd composed(3,6);
    composed.leftCols(3) = y.skew().evaluate();
    composed.rightCols(3) = w.skew().evaluate();

    assert(composed == horzcat.evaluate());

    assert(y.skew().getColumnDerivative(1, y).evaluate() == horzcat.getColumnDerivative(1,y).evaluate());
    assert(w.skew().getColumnDerivative(1, w).evaluate() == horzcat.getColumnDerivative(4,w).evaluate());

    levi::Expression vertcat = levi::Expression::Vertcat(y.skew(), w.skew(), "horzcat");

    Eigen::MatrixXd composedVert(6,3);
    composedVert.topRows(3) = y.skew().evaluate();
    composedVert.bottomRows(3) = w.skew().evaluate();

    assert(composedVert == vertcat.evaluate());

    assert(y.skew().getColumnDerivative(1, y).evaluate() == vertcat.getColumnDerivative(1,y).evaluate().topRows(3));
    assert(w.skew().getColumnDerivative(1, w).evaluate() == vertcat.getColumnDerivative(1,w).evaluate().bottomRows(3));


    return 0;

}
