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

    Eigen::MatrixXd A(Eigen::MatrixXd::Identity(3,3));
    Constant a(Eigen::MatrixXd::Identity(3,3), "a");
    Eigen::Matrix<double,3,3> B;
    B.setIdentity();
    Mutable b(3,3,"b");

    auto sum  = A + a;

    auto subtraction  = A - a;

    b = B;

    Expression c = a + b;

    auto squeeze = c.squeeze("squeezeTest");
    assert(squeeze.evaluate() == c.evaluate());

    c = c * a;
    std::cerr << c.name() << ": " <<std::endl << c.evaluate() <<std::endl;

    Variable x(3,"x"), y(3, "y");

    auto f = x + y;

    Eigen::Matrix<double,3,3> V;
    V.setIdentity();
    V *= 3;
    Constant v(V, "V");
    Expression g = v * c;
    Expression k = g * V;
    std::cerr << g.name() << " " <<std::endl << g.evaluate() <<std::endl;
    std::cerr << k.name() << " " <<std::endl << k.evaluate() <<std::endl;


    Scalar test(0.5);
    c = test* c;
    std::cerr << c.name() << ": " <<std::endl << c.evaluate() <<std::endl;

    c = 0.5* c + c * 0.5;
    std::cerr << c.name() << ": " <<std::endl << c.evaluate() <<std::endl;

    Expression d = 1.0 * c;

    assert(d.evaluate() == c.evaluate());

    d = 1.0 * d + c;

    assert(d.squeeze("squeezeTest2").evaluate() == d.evaluate());

    Eigen::VectorXd x_value(3);
    x_value.setRandom();
    x = x_value;
    ColumnExpression testX = g*x;
    std::cerr << testX.name() << " " <<std::endl << testX.evaluate() << std::endl <<"with x= " << x.evaluate()<<std::endl;
    x_value.setRandom();
    x = x_value;
    std::cerr << testX.name() << " " <<std::endl << testX.evaluate() << std::endl <<"with x= " << x.evaluate()<<std::endl;

    assert(x.evaluate() == x_value);

    auto row1 = c.row(0);

    row1.evaluate();

    assert(row1.evaluate() == c.evaluate().row(0));

    auto col1 = c.col(0);

    col1.evaluate();

    assert(col1.evaluate() == c.evaluate().col(0));

    auto testRow = test.row(0);

    auto testCol = test.col(0);

    auto cElement = c(0, 1);
    std::cerr << cElement.name() << " " <<std::endl << cElement.evaluate() <<std::endl;
    assert(cElement.evaluate() == c.evaluate()(0, 1));

    auto testElement = test(0,0);

    auto testBlock = c.block(1,1, 2, 2);
    std::cerr << testBlock.name() << " " <<std::endl << testBlock.evaluate() <<std::endl;
    assert(testBlock.evaluate() == c.evaluate().block(1,1,2,2));

    auto testSkew = x.skew();
    auto testVee = testSkew.vee();

    assert(x.evaluate() == testVee.evaluate());

    Variable z(1,"z");

    ScalarExpression zInverted = z.pow(-1);
    z = 2;
    assert(zInverted.evaluate() == 0.5);

    Expression zSquared = z/zInverted;

    assert(zSquared.evaluate()(0,0) == 4.0);

    assert(zSquared.getColumnDerivative(0, z).evaluate()(0,0) == 4.0);

    assert(zSquared.getColumnDerivative(0, z).getColumnDerivative(0, z).evaluate()(0,0) == 2.0);

    FixedSizeExpression<3,3> fixedExpr;
    fixedExpr = +v * Identity(1,1);
    assert(fixedExpr == v);

    Null zero(3,3);
    Expression testZero = zero;

    return 0;
}
