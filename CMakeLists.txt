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

find_package(iDynTree 0.11.103 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(levi REQUIRED)
find_package(matlogger2 REQUIRED)

set(CONSTRAINTS_HEADERS_DIR include/DynamicalPlannerPrivate/Constraints)

list(APPEND DPLANNER_PRIVATE_HEADERS ${CONSTRAINTS_HEADERS_DIR}/DynamicalConstraints.h
                                     ${CONSTRAINTS_HEADERS_DIR}/CoMPositionConstraint.h
                                     ${CONSTRAINTS_HEADERS_DIR}/QuaternionNormConstraint.h
                                     ${CONSTRAINTS_HEADERS_DIR}/ContactPositionConsistencyConstraint.h
                                     ${CONSTRAINTS_HEADERS_DIR}/FeetLateralDistanceConstraint.h
                                     ${CONSTRAINTS_HEADERS_DIR}/ContactFrictionConstraint.h
                                     ${CONSTRAINTS_HEADERS_DIR}/CentroidalMomentumConstraint.h
                                     ${CONSTRAINTS_HEADERS_DIR}/ContactForceControlConstraints.h
                                     ${CONSTRAINTS_HEADERS_DIR}/NormalVelocityControlConstraints.h
                                     ${CONSTRAINTS_HEADERS_DIR}/PlanarVelocityControlConstraints.h
                                     ${CONSTRAINTS_HEADERS_DIR}/DynamicalComplementarityConstraint.h
                                     include/DynamicalPlannerPrivate/Constraints.h)

set(COSTS_HEADERS_DIR include/DynamicalPlannerPrivate/Costs)

list(APPEND DPLANNER_PRIVATE_HEADERS ${COSTS_HEADERS_DIR}/FrameOrientationCost.h
                                     ${COSTS_HEADERS_DIR}/ForceMeanCost.h
                                     ${COSTS_HEADERS_DIR}/StaticTorquesCost.h
                                     ${COSTS_HEADERS_DIR}/SwingCost.h
                                     ${COSTS_HEADERS_DIR}/PhantomForcesCost.h
                                     ${COSTS_HEADERS_DIR}/MeanPointPositionCost.h
                                     ${COSTS_HEADERS_DIR}/CoMVelocityCost.h
                                     ${COSTS_HEADERS_DIR}/FootYawCost.h
                                     ${COSTS_HEADERS_DIR}/FeetDistanceCost.h
                                     include/DynamicalPlannerPrivate/Costs.h)

set(UTILITIES_DIR include/DynamicalPlannerPrivate/Utilities)

list(APPEND DPLANNER_PRIVATE_HEADERS ${UTILITIES_DIR}/VariablesLabeller.h
                                     ${UTILITIES_DIR}/QuaternionUtils.h
                                     ${UTILITIES_DIR}/SharedKinDynComputations.h
                                     ${UTILITIES_DIR}/CheckEqualVector.h
                                     ${UTILITIES_DIR}/HyperbolicSecant.h
                                     ${UTILITIES_DIR}/SmoothingFunction.h
                                     ${UTILITIES_DIR}/HyperbolicTangent.h
                                     ${UTILITIES_DIR}/TimelySharedKinDynComputations.h
                                     ${UTILITIES_DIR}/ExpressionsServer.h)

set(LEVI_UTILITIES_DIR include/DynamicalPlannerPrivate/Utilities/levi)

list(APPEND DPLANNER_PRIVATE_HEADERS ${LEVI_UTILITIES_DIR}/QuaternionExpressions.h
                                     ${LEVI_UTILITIES_DIR}/AdjointTransformExpression.h
                                     ${LEVI_UTILITIES_DIR}/RelativeJacobianExpression.h
                                     ${LEVI_UTILITIES_DIR}/RelativeQuaternionExpression.h
                                     ${LEVI_UTILITIES_DIR}/RelativePositionExpression.h
                                     ${LEVI_UTILITIES_DIR}/TransformExpression.h
                                     ${LEVI_UTILITIES_DIR}/CoMInBaseExpression.h
                                     ${LEVI_UTILITIES_DIR}/RelativeVelocityExpression.h
                                     ${LEVI_UTILITIES_DIR}/QuaternionErrorExpression.h
                                     ${LEVI_UTILITIES_DIR}/AbsoluteVelocityExpression.h
                                     ${LEVI_UTILITIES_DIR}/MomentumInBaseExpression.h)

set(DPLANNER_PRIVATE_SOURCES src/private/VariablesLabeller.cpp
                             src/private/DynamicalConstraints.cpp
                             src/private/QuaternionUtils.cpp
                             src/private/CoMPositionConstraint.cpp
                             src/private/QuaternionNormConstraint.cpp
                             src/private/ContactPositionConsistencyConstraint.cpp
                             src/private/SharedKinDynComputations.cpp
                             src/private/FeetLateralDistanceConstraint.cpp
                             src/private/ContactFrictionConstraint.cpp
                             src/private/CentroidalMomentumConstraint.cpp
                             src/private/HyperbolicSecant.cpp
                             src/private/ContactForceControlConstraints.cpp
                             src/private/NormalVelocityControlConstraints.cpp
                             src/private/PlanarVelocityControlConstraints.cpp
                             src/private/FrameOrientationCost.cpp
                             src/private/ForceMeanCost.cpp
                             src/private/StaticTorquesCost.cpp
                             src/private/SmoothingFunction.cpp
                             src/private/HyperbolicTangent.cpp
                             src/private/DynamicalComplementarityConstraint.cpp
                             src/private/TimelySharedKinDynComputations.cpp
                             src/private/SwingCost.cpp
                             src/private/PhantomForcesCost.cpp
                             src/private/MeanPointPositionCost.cpp
                             src/private/CoMVelocityCost.cpp
                             src/private/QuaternionExpressions.cpp
                             src/private/AdjointTransformExpression.cpp
                             src/private/RelativeJacobianExpression.cpp
                             src/private/RelativeQuaternionExpression.cpp
                             src/private/RelativePositionExpression.cpp
                             src/private/TransformExpression.cpp
                             src/private/CoMInBaseExpression.cpp
                             src/private/RelativeVelocityExpression.cpp
                             src/private/ExpressionsServer.cpp
                             src/private/QuaternionErrorExpression.cpp
                             src/private/AbsoluteVelocityExpression.cpp
                             src/private/MomentumInBaseExpression.cpp
                             src/private/FootYawCost.cpp
                             src/private/FeetDistanceCost.cpp)


add_library(DynamicalPlannerPrivate ${DPLANNER_PRIVATE_HEADERS} ${DPLANNER_PRIVATE_SOURCES})
target_include_directories(DynamicalPlannerPrivate PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>)
target_include_directories(DynamicalPlannerPrivate PRIVATE ${EIGEN3_INCLUDE_DIR})
target_link_libraries(DynamicalPlannerPrivate PUBLIC ${iDynTree_LIBRARIES})
target_link_libraries(DynamicalPlannerPrivate PUBLIC levi::levi)

set(DPLANNER_HEADERS include/DynamicalPlanner/Settings.h
                     include/DynamicalPlanner/Solver.h
                     include/DynamicalPlanner/State.h
                     include/DynamicalPlanner/Control.h
                     include/DynamicalPlanner/Visualizer.h
                     include/DynamicalPlanner/RectangularFoot.h
                     include/DynamicalPlanner/Logger.h)

set(DPLANNER_SOURCES src/Settings.cpp
                     src/Solver.cpp
                     src/State.cpp
                     src/Control.cpp
                     src/RectangularFoot.cpp
                     src/Visualizer.cpp
                     src/Logger.cpp)

add_library(DynamicalPlanner ${DPLANNER_HEADERS} ${DPLANNER_SOURCES})
target_include_directories(DynamicalPlanner PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories(DynamicalPlanner PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(DynamicalPlannerPrivate PRIVATE ${EIGEN3_INCLUDE_DIR})
target_link_libraries(DynamicalPlanner PRIVATE DynamicalPlannerPrivate)
target_link_libraries(DynamicalPlanner PRIVATE matlogger2::matlogger2)
target_link_libraries(DynamicalPlanner PUBLIC ${iDynTree_LIBRARIES})

include(AddUninstallTarget)

option(BUILD_TESTING "Create tests using CMake" OFF)

option(RUN_VALGRIND_TESTS "Run tests with Valgrind" FALSE)
mark_as_advanced(RUN_VALGRIND_TESTS)
if(BUILD_TESTING)
   include( CTest )
   enable_testing()
endif()


if(BUILD_TESTING)
    add_subdirectory(test)
endif()
