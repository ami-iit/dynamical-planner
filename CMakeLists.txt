# Copyright (C) 2013-2018 Fondazione Istituto Italiano di Tecnologia
#
# Licensed under either the GNU Lesser General Public License v3.0 :
# https://www.gnu.org/licenses/lgpl-3.0.html
# or the GNU Lesser General Public License v2.1 :
# https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
# at your option.

cmake_minimum_required(VERSION 3.5)

project(DynamicalPlanner CXX)

include(GNUInstallDirs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(MSVC)
    set(CMAKE_DEBUG_POSTFIX "d")
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)

# Enable RPATH support for installed binaries and libraries
include(AddInstallRPATHSupport)
add_install_rpath_support(BIN_DIRS "${CMAKE_INSTALL_FULL_LIBDIR}"
                          LIB_DIRS "${CMAKE_INSTALL_FULL_BINDIR}"
                          INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}"
                          USE_LINK_PATH)

# Encourage user to specify a build type (e.g. Release, Debug, etc.), otherwise set it to Release.
if(NOT CMAKE_CONFIGURATION_TYPES)
    if(NOT CMAKE_BUILD_TYPE)
        message(STATUS "Setting build type to 'Release' as none was specified.")
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release")
    endif()
endif()


find_package(iDynTree 0.10.0 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})


set(DPLANNER_HEADERS )

set(DPLANNER_PRIVATE_HEADERS include/DynamicalPlannerPrivate/VariablesLabeller.h
                             include/DynamicalPlannerPrivate/DynamicalConstraints.h
                             include/DynamicalPlannerPrivate/QuaternionUtils.h
                             include/DynamicalPlannerPrivate/CoMPositionConstraint.h
                             include/DynamicalPlannerPrivate/QuaternionNormConstraint.h
                             include/DynamicalPlannerPrivate/ContactPositionConsistencyConstraint.h
                             include/DynamicalPlannerPrivate/SharedKinDynComputations.h
                             include/DynamicalPlannerPrivate/CheckEqualVector.h
                             include/DynamicalPlannerPrivate/FeetLateralDistanceConstraint.h
                             include/DynamicalPlannerPrivate/ContactFrictionConstraint.h
                             include/DynamicalPlannerPrivate/CentroidalMomentumConstraint.h
                             include/DynamicalPlannerPrivate/HyperbolicSecant.h
                             include/DynamicalPlannerPrivate/ContactForceControlConstraints.h
                             include/DynamicalPlannerPrivate/ContactVelocityControlConstraints.h
                             include/DynamicalPlannerPrivate/Constraints.h
                             include/DynamicalPlannerPrivate/FrameOrientationCost.h
                             include/DynamicalPlannerPrivate/ForceMeanCost.h
                             include/DynamicalPlannerPrivate/Costs.h)

set(DPLANNER_SOURCES src/VariablesLabeller.cpp
                     src/DynamicalConstraints.cpp
                     src/QuaternionUtils.cpp
                     src/CoMPositionConstraint.cpp
                     src/QuaternionNormConstraint.cpp
                     src/ContactPositionConsistencyConstraint.cpp
                     src/SharedKinDynComputations.cpp
                     src/FeetLateralDistanceConstraint.cpp
                     src/ContactFrictionConstraint.cpp
                     src/CentroidalMomentumConstraint.cpp
                     src/HyperbolicSecant.cpp
                     src/ContactForceControlConstraints.cpp
                     src/ContactVelocityControlConstraints.cpp
                     src/FrameOrientationCost.cpp
                     src/ForceMeanCost.cpp)

add_library(DynamicalPlanner ${DPLANNER_HEADERS} ${DPLANNER_PRIVATE_HEADERS} ${DPLANNER_SOURCES})
target_include_directories(DynamicalPlanner PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories(DynamicalPlanner PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(DynamicalPlanner ${iDynTree_LIBRARIES})


include(AddUninstallTarget)

option(BUILD_TESTING "Create tests using CMake" OFF)

option(RUN_VALGRIND_TESTS "Run tests with Valgrind" FALSE)
mark_as_advanced(RUN_VALGRIND_TESTS)
if(BUILD_TESTING)
   include( CTest )
   enable_testing()

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
endif()


if(BUILD_TESTING)
    add_subdirectory(test)
endif()
