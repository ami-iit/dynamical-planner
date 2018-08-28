/*
 * Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
 * Authors: Stefano Dafarra
 * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 *
 */

#include <DynamicalPlannerPrivate/VariablesLabeller.h>
#include <iDynTree/Core/TestUtils.h>
#include <iDynTree/Core/VectorDynSize.h>
#include <iDynTree/Core/EigenHelpers.h>

int main()
{
    using namespace DynamicalPlanner::Private;
    VariablesLabeller variables;

    iDynTree::VectorDynSize fullVector(9), testVector(9), testPart(3);
    iDynTree::getRandomVector(fullVector);

    ASSERT_IS_TRUE(variables.addLabel("part1", 3));
    ASSERT_IS_TRUE(variables.addLabel("part2", 3));
    ASSERT_IS_TRUE(variables.addLabel("part3", 3));

    variables = fullVector;
    testVector = variables.values();
    ASSERT_IS_TRUE(variables.size() == 9);
    ASSERT_EQUAL_VECTOR(testVector, fullVector);
    iDynTree::IndexRange range;

    range = variables.getIndexRange("part2");
    ASSERT_IS_TRUE(range.isValid());
    testPart = variables("part2");
    testVector.resize(3);
    iDynTree::toEigen(testVector) = iDynTree::toEigen(fullVector).segment(range.offset, range.size);
    ASSERT_EQUAL_VECTOR(testVector, testPart);
    iDynTree::toEigen(variables("part2")).setZero();

    range = variables.getIndexRange("part1");
    ASSERT_IS_TRUE(range.isValid());
    testPart = variables("part1");
    testVector.resize(3);
    iDynTree::toEigen(testVector) = iDynTree::toEigen(fullVector).segment(range.offset, range.size);
    ASSERT_EQUAL_VECTOR(testVector, testPart);
    iDynTree::toEigen(variables("part1")).setZero();


    range = variables.getIndexRange("part3");
    ASSERT_IS_TRUE(range.isValid());
    testPart = variables(range);
    testVector.resize(3);
    iDynTree::toEigen(testVector) = iDynTree::toEigen(fullVector).segment(range.offset, range.size);
    ASSERT_EQUAL_VECTOR(testVector, testPart);
    iDynTree::toEigen(variables("part3")).setZero();

    iDynTree::toEigen(fullVector).setZero();
    testVector.resize(9);
    testVector = variables.values();
    ASSERT_EQUAL_VECTOR(testVector, fullVector);

    ASSERT_IS_TRUE(variables.listOfLabels().size() == 3);


    return EXIT_SUCCESS;
}
