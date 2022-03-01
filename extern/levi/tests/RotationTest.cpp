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

    FixedSizeVariable<4> normalizedQuaternion = (quaternion/(quaternion.transpose() * quaternion).pow(0.5)).asVariable();

    FixedSizeExpression<3,3> skewQuaternion = normalizedQuaternion.block<3,1>(1,0).skew();

    FixedSizeExpression<3,3> twoSkewQuaternion = 2.0 * skewQuaternion;

    FixedSizeVariable<3> x(3, "x");

    FixedSizeExpression<3,3> rotation;

    //Rodriguez formula
    rotation = FixedSizeIdentity<3,3>() + normalizedQuaternion(0,0) * twoSkewQuaternion + twoSkewQuaternion * skewQuaternion;

    std::cerr << rotation.name() << std::endl;

    Eigen::Vector4d quaternionValue;

    quaternionValue.setZero();
    quaternionValue[0] = 1.0;

    quaternion = quaternionValue;

    Eigen::MatrixXd testRotation(3,3);
    testRotation.setIdentity();

    assert(rotation.evaluate() == testRotation);

    quaternionValue = Eigen::Vector4d::Random();

    quaternion = quaternionValue;
    assert(rotation.isNew());

    Eigen::Vector4d quaternionValueNormalized = quaternionValue.normalized();

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> skew;
    skew << 0.0, -quaternionValueNormalized[3], quaternionValueNormalized[2],
          quaternionValueNormalized[3], 0.0, -quaternionValueNormalized[1],
         -quaternionValueNormalized[2], quaternionValueNormalized[1], 0.0;

    testRotation = Eigen::MatrixXd::Identity(3,3) + 2.0 * quaternionValueNormalized(0) * skew + 2.0 * skew * skew;

    rotation.evaluate();
    assert(rotation.evaluate() == testRotation);

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

    begin = std::chrono::steady_clock::now();
    testSpeed = squeezed.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (squeezed): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    assert(squeezed.evaluate() == rotation.evaluate());


    //-------------------------Validation of first derivative

    double perturbationValue = 1e-3;
    Eigen::Vector4d quaternionPerturbed, perturbation;
    Eigen::Vector3d vector, output, perturbedOutput, firstOrderTaylor;
    vector = Eigen::Vector3d::Random() * 10.0;
    Eigen::Matrix<double, 3, 3> perturbedRotation;

    Eigen::Vector4d maxQuaternion;
    maxQuaternion.setConstant(1.0);
    Eigen::Vector4d minQuaternion;
    minQuaternion.setConstant(-1);
    minQuaternion(0) = 0;

    Expression rotatedVector = rotation * x;

    x = vector;

    output = rotatedVector.evaluate();

    begin = std::chrono::steady_clock::now();
    Expression rotatedVectorDerivative = rotatedVector.getColumnDerivative(0, normalizedQuaternion) * normalizedQuaternion.getColumnDerivative(0, quaternion);
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (compute first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    Eigen::MatrixXd derivativeValue(rotatedVectorDerivative.rows(), rotatedVectorDerivative.cols());
    begin = std::chrono::steady_clock::now();
    derivativeValue = rotatedVectorDerivative.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (evaluate first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

    auto squeezedDerivative = rotatedVectorDerivative.squeeze("squeezedDerivative");
    Eigen::MatrixXd squeezedDerivativeValue(rotatedVectorDerivative.rows(), rotatedVectorDerivative.cols());

    begin = std::chrono::steady_clock::now();
    squeezedDerivativeValue = squeezedDerivative.evaluate();
    end= std::chrono::steady_clock::now();
    std::cout << "Elapsed time ms (squeezed first derivative): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;
    assert((squeezedDerivativeValue - derivativeValue).cwiseAbs().maxCoeff() < 1e-10);


//    std::cerr << rotatedVectorDerivative.name() << std::endl;


    // Test separetly the derivative of quaternion
    for (unsigned int i = 0; i < 4; i++)
    {
        quaternionPerturbed = quaternionValue;
        quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

        //ensure validity of obtained quaternion even if the quaternion is no more a step different than
        //the original
//        quaternionPerturbed = quaternionPerturbed.array().min(maxQuaternion.array());
//        quaternionPerturbed = quaternionPerturbed.array().max(minQuaternion.array());
//        quaternionPerturbed.normalize();

        perturbation = quaternionPerturbed - quaternionValue;

        quaternion = quaternionPerturbed;

        perturbedOutput = rotatedVector.evaluate();

        firstOrderTaylor = derivativeValue * perturbation;

        firstOrderTaylor += output;

        assert ((firstOrderTaylor - perturbedOutput).cwiseAbs().maxCoeff() < perturbationValue/50.0);

    }

    //-------------------------Validation of second derivative

    for (Eigen::Index j = 0 ; j < rotatedVectorDerivative.cols(); ++j) {
        quaternion = quaternionValue;

        output = rotatedVectorDerivative.col(j).evaluate();

        begin = std::chrono::steady_clock::now();
        Expression rotatedVectorDoubleDerivative = rotatedVectorDerivative.getColumnDerivative(j, quaternion);
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (compute second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

        auto squeezedSecondDerivative = rotatedVectorDoubleDerivative.squeeze("doubleSqueezedDerivative" + std::to_string(j));

        Eigen::MatrixXd squeezedOutput(rotatedVectorDoubleDerivative.rows(), rotatedVectorDoubleDerivative.cols());

        begin = std::chrono::steady_clock::now();
        squeezedOutput = squeezedSecondDerivative.evaluate();
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (evaluate (squeeze) second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;


        Eigen::MatrixXd doubleDerivativeValue(rotatedVectorDoubleDerivative.rows(), rotatedVectorDoubleDerivative.cols());

//        std::cerr << rotatedVectorDoubleDerivative.name() << std::endl;

        quaternion = quaternionValue;

        begin = std::chrono::steady_clock::now();
        doubleDerivativeValue = rotatedVectorDoubleDerivative.evaluate();
        end= std::chrono::steady_clock::now();
        std::cout << "Elapsed time ms (evaluate second derivative, column " << j<<"): " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000.0) <<std::endl;

        // Test separetly the derivative of quaternion
        for (unsigned int i = 0; i < 4; i++)
        {
            quaternionPerturbed = quaternionValue;
            quaternionPerturbed(i) = quaternionPerturbed(i) + perturbationValue;

            //ensure validity of obtained quaternion even if the quaternion is no more a step different than
            //the original
//            quaternionPerturbed = quaternionPerturbed.array().min(maxQuaternion.array());
//            quaternionPerturbed = quaternionPerturbed.array().max(minQuaternion.array());
    //        quaternionPerturbed.normalize();

            perturbation = quaternionPerturbed - quaternionValue;

            quaternion = quaternionPerturbed;

            perturbedOutput = rotatedVectorDerivative.evaluate().col(j);

            firstOrderTaylor = doubleDerivativeValue * perturbation;

            firstOrderTaylor += output;

            assert ((firstOrderTaylor - perturbedOutput).cwiseAbs().maxCoeff() < perturbationValue/20.0);

        }
    }

    return 0;
}
