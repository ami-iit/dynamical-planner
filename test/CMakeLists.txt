# Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
#
# Licensed under either the GNU Lesser General Public License v3.0 :
# https://www.gnu.org/licenses/lgpl-3.0.html
# or the GNU Lesser General Public License v2.1 :
# https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
# at your option.

include_directories(${CMAKE_CURRENT_BINARY_DIR}/data)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/URDFdir.h.in" "${CMAKE_CURRENT_BINARY_DIR}/data/URDFdir.h" @ONLY)

macro(add_dp_test classname)
    set(testsrc ${classname}Test.cpp)
    set(testbinary ${classname}UnitTest)
    set(testname   UnitTest${classname})
    add_executable(${testbinary} ${testsrc})
    target_include_directories(${testbinary} PRIVATE ${IDYNTREE_TREE_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(${testbinary} PRIVATE DynamicalPlanner)
    add_test(NAME ${testname} COMMAND ${testbinary})

    if(RUN_VALGRIND_TESTS)
        add_test(NAME memcheck_${testname} COMMAND ${MEMCHECK_COMMAND_COMPLETE} $<TARGET_FILE:${testbinary}>)
    endif()
endmacro()

add_dp_test(VariablesLabeller)
add_dp_test(QuaternionDerivative)
add_dp_test(ConstraintsDerivative)
add_dp_test(MomentumDerivative)

