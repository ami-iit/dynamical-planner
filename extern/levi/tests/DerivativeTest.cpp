/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <levi/levi.h>
#include <iostream>

int main() {
    using namespace Eigen;
    using namespace levi;

    Variable x(3, "x"), y(3, "y");

    Expression xDerivative = x.getColumnDerivative(0, x);

    std::cerr << xDerivative.name() << std::endl <<xDerivative.evaluate() << std::endl;

    Expression yDerivative = y.getColumnDerivative(0, x);

    std::cerr << yDerivative.name() << std::endl << yDerivative.evaluate() << std::endl;

    Expression variableSum = x + y;

    auto sumDerivative = variableSum.getColumnDerivative(0, x);

    std::cerr << sumDerivative.name() << std::endl << sumDerivative.evaluate() << std::endl;

    MatrixXd A(MatrixXd::Identity(3,3));
    Constant a(MatrixXd::Identity(3,3), "a");
    Matrix<double,3,3> B;
    B.setIdentity();
    Mutable b(3,3,"b");

    b = B;

    Expression expr = 2.0 * a * b * x + b*y;

    Expression exprDerivative = expr.getColumnDerivative(0, x);

    std::cerr << exprDerivative.name() << std::endl << exprDerivative.evaluate() << std::endl;

    return 0;
}
