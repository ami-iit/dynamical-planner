# SPDX-FileCopyrightText: Fondazione Istituto Italiano di Tecnologia (IIT)
# SPDX-License-Identifier: BSD-3-Clause

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Enable C++14
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
    set(CMAKE_CXX_STANDARD 17) #Workaround to fix Eigen alignement problems when using fixed size matrices
else()
    set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Turn on testing.
option(RUN_VALGRIND_TESTS "Run tests with Valgrind" FALSE)
mark_as_advanced(RUN_VALGRIND_TESTS)

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

macro(add_levi_test classname)
    set(testsrc ${classname}Test.cpp)
    set(testbinary ${classname}UnitTest)
    set(testname   UnitTest${classname})
    add_executable(${testbinary} ${testsrc})
    target_include_directories(${testbinary} PRIVATE ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(${testbinary} PRIVATE levi::levi)
    add_test(NAME ${testname} COMMAND ${testbinary})

    if(RUN_VALGRIND_TESTS)
        add_test(NAME memcheck_${testname} COMMAND ${MEMCHECK_COMMAND_COMPLETE} $<TARGET_FILE:${testbinary}>)
    endif()
endmacro()

add_levi_test(Composition)
add_levi_test(Derivative)
add_levi_test(AdvancedConstructors)
add_levi_test(Cache)
add_levi_test(Dependencies)
add_levi_test(Rotation)

option(ENABLE_COMPILED_ROTATION_TEST "Enable the test on the automatic code generation and compilation" OFF)

if (ENABLE_COMPILED_ROTATION_TEST)
    add_levi_test(CompiledRotation)
endif()

