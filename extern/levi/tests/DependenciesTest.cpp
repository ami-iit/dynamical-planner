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

    Variable x(3, "x");
    Variable y(4, "y");

    auto dependecyList = make_dependencies_set(x,y, Null(3,3));

    size_t size = dependecyList.size();
    assert(size == 3);

    bool isFirstNew = dependecyList.isNew(levi::size_t_value<0>());
    assert(isFirstNew);

    bool isNew = dependecyList.areNew();

    assert(isNew);

    Expression test = dependecyList.get(levi::size_t_value<0>());

    dependecyList.evaluate(levi::size_t_value<0>());

    isFirstNew = dependecyList.isNew(levi::size_t_value<0>());

    assert(!isFirstNew);

    bool isNull = dependecyList.get(levi::size_t_value<2>()).isNull();
    assert(isNull);

    Variable z(3, "z");
    Variable k(4, "k");

    auto otherDependecyList = make_dependencies_set(z,k, Null(3,3));
    otherDependecyList.evaluate(levi::size_t_value<1>());

    dependecyList = otherDependecyList;

    return 0;
}
