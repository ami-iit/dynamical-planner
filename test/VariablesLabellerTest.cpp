/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <private/VariablesLabeller.h>
#include <private/SpanUtils.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/EigenHelpers.h>

int main()
{
    VariablesLabeller variables;

    iDynTree::VectorDynSize fullVector(9), testVector(9), testPart(3);
    iDynTree::getRandomVector(fullVector);

    ASSERT_IS_TRUE(variables.addLabel("part1", 3));
    ASSERT_IS_TRUE(variables.addLabel("part2", 3));
    ASSERT_IS_TRUE(variables.addLabel("part3", 3));

    spanToEigen(variables.values()) = iDynTree::toEigen(fullVector);
    iDynTree::toEigen(testVector) = spanToEigen(variables.values());
    ASSERT_EQUAL_VECTOR(testVector, fullVector);
    iDynTree::IndexRange range;

    range = variables.getIndexRange("part2");
    ASSERT_IS_TRUE(range.isValid());
    iDynTree::toEigen(testPart) = spanToEigen(variables.values("part2"));
    testVector.resize(3);
    iDynTree::toEigen(testVector) = iDynTree::toEigen(fullVector).segment(range.offset, range.size);
    ASSERT_EQUAL_VECTOR(testVector, testPart);
    spanToEigen(variables.values("part2")).setZero();

    range = variables.getIndexRange("part1");
    ASSERT_IS_TRUE(range.isValid());
    iDynTree::toEigen(testPart) = spanToEigen(variables.values("part1"));
    testVector.resize(3);
    iDynTree::toEigen(testVector) = iDynTree::toEigen(fullVector).segment(range.offset, range.size);
    ASSERT_EQUAL_VECTOR(testVector, testPart);
    spanToEigen(variables.values("part1")).setZero();


    range = variables.getIndexRange("part3");
    ASSERT_IS_TRUE(range.isValid());
    iDynTree::toEigen(testPart) = spanToEigen(variables.values("part3"));
    testVector.resize(3);
    iDynTree::toEigen(testVector) = iDynTree::toEigen(fullVector).segment(range.offset, range.size);
    ASSERT_EQUAL_VECTOR(testVector, testPart);
    spanToEigen(variables.values("part3")).setZero();

    iDynTree::toEigen(fullVector).setZero();
    testVector.resize(9);
    iDynTree::toEigen(testVector) = spanToEigen(variables.values());
    ASSERT_EQUAL_VECTOR(testVector, fullVector);


    return EXIT_SUCCESS;
}
