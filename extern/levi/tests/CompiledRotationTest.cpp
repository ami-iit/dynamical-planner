/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */


#include <levi/levi.h>
#include <Eigen/Geometry>
#include <iostream>
#include <chrono>


int main() {
    using namespace levi;

    FixedSizeVariable<4> quaternion(4, "q");

    FixedSizeExpression<4, 1> normalizedQuaternion = (quaternion/(quaternion.transpose() * quaternion).pow(0.5));

    FixedSizeExpression<3,3> skewQuaternion = normalizedQuaternion.block<3,1>(1,0).skew();

    FixedSizeExpression<3,3> twoSkewQuaternion = 2.0 * skewQuaternion;

    FixedSizeVariable<3> x(3, "x");

    FixedSizeExpression<3,3> rotation;

    //Rodriguez formula
    rotation = FixedSizeIdentity<3,3>() + normalizedQuaternion(0,0) * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;

    std::cerr << rotation.name() << std::endl;

    Eigen::Vector4d quaternionValue;

    quaternionValue = Eigen::Vector4d::Random();

    quaternion = quaternionValue;
    assert(rotation.isNew());

    Eigen::Vector4d quaternionValueNormalized = quaternionValue.normalized();

    Eigen::Quaterniond quaternionEigen;
    quaternionEigen.w() = quaternionValueNormalized[0];
    quaternionEigen.x() = quaternionValueNormalized[1];
    quaternionEigen.y() = quaternionValueNormalized[2];
    quaternionEigen.z() = quaternionValueNormalized[3];

    assert(rotation.evaluate() == quaternionEigen.toRotationMatrix());

    Eigen::Matrix<double, 3, 3> testSpeed;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    testSpeed = quaternionEigen.toRotationMatrix();
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (eigen): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    quaternion = quaternionValue;

    assert(rotation.isNew());
    begin = std::chrono::steady_clock::now();
    testSpeed = rotation.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (evaluate): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    auto squeezed = rotation.squeeze("SqueezedRotation");

    quaternion = quaternionValue;

    begin = std::chrono::steady_clock::now();
    testSpeed = squeezed.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (squeezed): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    assert(squeezed.evaluate() == rotation.evaluate());


    auto compiled = rotation.compile("SqueezedRotation");

    quaternion = quaternionValue;

    begin = std::chrono::steady_clock::now();
    testSpeed = compiled.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (compiled): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    assert(compiled.evaluate() == rotation.evaluate());


    //-------------------------Validation of first derivative

    Eigen::Vector3d vector;
    vector = Eigen::Vector3d::Random() * 10.0;

    Expression rotatedVector = rotation * x;

    quaternion = quaternionValue;
    x = vector;

    begin = std::chrono::steady_clock::now();
    Expression rotatedVectorDerivative = rotatedVector.getColumnDerivative(0, quaternion);
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (compute first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    quaternion = quaternionValue;
    x = vector;

    Eigen::MatrixXd derivativeValue(rotatedVectorDerivative.rows(), rotatedVectorDerivative.cols());
    begin = std::chrono::steady_clock::now();
    derivativeValue = rotatedVectorDerivative.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (evaluate first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    auto squeezedDerivative = rotatedVectorDerivative.squeeze("squeezedDerivative");
    Eigen::MatrixXd squeezedDerivativeValue(rotatedVectorDerivative.rows(), rotatedVectorDerivative.cols());

    quaternion = quaternionValue;
    x = vector;

    begin = std::chrono::steady_clock::now();
    squeezedDerivativeValue = squeezedDerivative.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (squeezed first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
    assert((squeezedDerivativeValue - derivativeValue).cwiseAbs().maxCoeff() < 1e-10);

    auto compiledDerivative = rotatedVectorDerivative.compile("squeezedDerivative");
    Eigen::MatrixXd compiledDerivativeValue(rotatedVectorDerivative.rows(), rotatedVectorDerivative.cols());

    x = vector;
    quaternion = quaternionValue;

    begin = std::chrono::steady_clock::now();
    compiledDerivativeValue = compiledDerivative.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (compiled first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
    assert((compiledDerivativeValue - derivativeValue).cwiseAbs().maxCoeff() < 1e-10);


    //-------------------------Validation of second derivative

    MultipleExpressionsMap<Eigen::Matrix<double, 3, 4>> expressions;
    for (Eigen::Index j = 0 ; j < rotatedVectorDerivative.cols(); ++j) {

        x = vector;
        quaternion = quaternionValue;

        begin = std::chrono::steady_clock::now();
        Expression rotatedVectorDoubleDerivative = rotatedVectorDerivative.getColumnDerivative(j, quaternion);
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (compute second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

        Eigen::MatrixXd doubleDerivativeValue(rotatedVectorDoubleDerivative.rows(), rotatedVectorDoubleDerivative.cols());

        //        std::cerr << rotatedVectorDoubleDerivative.name() << std::endl;

        quaternion = quaternionValue;
        x = vector;

        begin = std::chrono::steady_clock::now();
        doubleDerivativeValue = rotatedVectorDoubleDerivative.evaluate();
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (evaluate second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

        auto squeezedSecondDerivative = rotatedVectorDoubleDerivative.squeeze("doubleSqueezedDerivative" + std::to_string(j));

        Eigen::MatrixXd squeezedOutput(rotatedVectorDoubleDerivative.rows(), rotatedVectorDoubleDerivative.cols());

        quaternion = quaternionValue;
        x = vector;

        begin = std::chrono::steady_clock::now();
        squeezedOutput = squeezedSecondDerivative.evaluate();
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (evaluate (squeeze) second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
        assert((squeezedOutput - doubleDerivativeValue).cwiseAbs().maxCoeff() < 1e-10);



        auto compiledSecondDerivative = rotatedVectorDoubleDerivative.compile("doubleSqueezedDerivative" + std::to_string(j));

        Eigen::MatrixXd compiledOutput(rotatedVectorDoubleDerivative.rows(), rotatedVectorDoubleDerivative.cols());

        quaternion = quaternionValue;
        x = vector;

        begin = std::chrono::steady_clock::now();
        compiledOutput = compiledSecondDerivative.evaluate();
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (evaluate (compile) second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
        assert((compiledOutput - doubleDerivativeValue).cwiseAbs().maxCoeff() < 1e-10);


        expressions["Col" + std::to_string(j)] = rotatedVectorDoubleDerivative;
    }

    auto multiExpr = levi::CompileMultipleExpressions(expressions, "RotatedVectorHessian");
    auto output = multiExpr->evaluate();

    auto multiSqueezed = levi::SqueezeMultipleExpressions(expressions);
    auto outputMultiSqueeze = multiSqueezed->evaluate();

    x = vector;
    quaternion = quaternionValue;

    begin = std::chrono::steady_clock::now();
    output = multiExpr->evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (evaluate (compile) second derivative, multiple): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    x = vector;
    quaternion = quaternionValue;

    begin = std::chrono::steady_clock::now();
    outputMultiSqueeze = multiSqueezed->evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (evaluate (squeezed) second derivative, multiple): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    Eigen::MatrixXd doubleDerivativeValue(3, 4);

    x = vector;
    quaternion = quaternionValue;

    begin = std::chrono::steady_clock::now();
    for (auto& expression: expressions) {
        doubleDerivativeValue = expression.second.evaluate();
    }
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (evaluate second derivatives, serial): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    for (auto& expression: expressions) {
        doubleDerivativeValue = expression.second.evaluate();
        assert((output[expression.first] - doubleDerivativeValue).cwiseAbs().maxCoeff() < 1e-10);
        assert((outputMultiSqueeze[expression.first] - doubleDerivativeValue).cwiseAbs().maxCoeff() < 1e-10);
    }

    return 0;
    }
