# Copyright (C) 2018 Fondazione Istituto Italiano di Tecnologia
#
# Licensed under either the GNU Lesser General Public License v3.0 :
# https://www.gnu.org/licenses/lgpl-3.0.html
# or the GNU Lesser General Public License v2.1 :
# https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
# at your option.

include_directories(${CMAKE_CURRENT_BINARY_DIR}/data)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/URDFdir.h.in" "${CMAKE_CURRENT_BINARY_DIR}/data/URDFdir.h" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/FolderPath.h.in" "${CMAKE_CURRENT_BINARY_DIR}/FolderPath.h" @ONLY)
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/SavedVideos")

# adding support for checking the tests with valgrind
if(RUN_VALGRIND_TESTS)
     find_package(Valgrind REQUIRED)
     if(VALGRIND_FOUND)
         set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_PROGRAM})
         set(MEMORYCHECK_COMMAND ${VALGRIND_PROGRAM})
         if (APPLE)
             set(MEMORYCHECK_SUPPRESSIONS "--suppressions=${PROJECT_SOURCE_DIR}/cmake/valgrind-macos.supp")
         else ()
             set(MEMORYCHECK_SUPPRESSIONS "")
         endif ()
         set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1 ${MEMORYCHECK_SUPPRESSIONS}"  CACHE STRING "Options to pass to the memory checker")
         mark_as_advanced(MEMORYCHECK_COMMAND_OPTIONS)
         set(MEMCHECK_COMMAND_COMPLETE "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
         separate_arguments(MEMCHECK_COMMAND_COMPLETE)
     endif()
endif()

macro(add_dp_test classname)
    set(testsrc ${classname}Test.cpp)
    set(testbinary ${classname}UnitTest)
    set(testname   UnitTest${classname})
    add_executable(${testbinary} ${testsrc})
    target_include_directories(${testbinary} PRIVATE {EIGEN3_INCLUDE_DIR})
    target_link_libraries(${testbinary} PRIVATE DynamicalPlanner)
    target_link_libraries(${testbinary} PRIVATE DynamicalPlannerPrivate)
    add_test(NAME ${testname} COMMAND ${testbinary})

    if(RUN_VALGRIND_TESTS)
        add_test(NAME memcheck_${testname} COMMAND ${MEMCHECK_COMMAND_COMPLETE} $<TARGET_FILE:${testbinary}>)
    endif()
endmacro()

add_dp_test(VariablesLabeller)
add_dp_test(QuaternionDerivative)
add_dp_test(ConstraintsDerivative)
add_dp_test(MomentumDerivative)
add_dp_test(CostsDerivative)
add_dp_test(StaticForcesDerivative)
add_dp_test(BoundsAndSparsity)
add_dp_test(RectangularFoot)
add_dp_test(Solver)
add_dp_test(SolverForComparisons)
add_dp_test(Performance)
add_dp_test(leviExpressions)
add_dp_test(Transcription)
add_dp_test(Logger)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/meshes" DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
