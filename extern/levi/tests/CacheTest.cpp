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

    Variable x(2, "x");

    Expression c = x + x.col(0);
    Expression c1 = c;

    Eigen::Vector2d x_value;

    for (size_t i = 0; i < 10; ++i) {
        x_value = Eigen::Vector2d::Random();
        x = x_value;

        assert(c.isNew() && c1.isNew());
        assert(c.evaluate() == (2.0 * x_value));

        assert(!c.isNew());
        assert(c1.isNew());
        assert(c1.evaluate() == (2.0 * x_value));
        assert(!c1.isNew());
    }

    for (size_t i = 0; i < 1000; ++i) {
        Expression newExpression; //test that the register does not grow indefinitely (this can be actually checked only with gdb)

        newExpression = c;

        assert(newExpression.isNew());
    }

    return 0;
}

