cmake_minimum_required(VERSION 3.5)

# ouptut paths
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# Build shared libs
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

if(MSVC)
    set(CMAKE_DEBUG_POSTFIX "d")
endif()

option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)


# Enable C++14
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# include macros for warnings
include(AddWarningsConfigurationToTargets)
include(CMakePackageConfigHelpers)

project(zupply
        LANGUAGES CXX)

include(GNUInstallDirs)

option(ENABLE_RPATH "Enable RPATH for this library" ON)
mark_as_advanced(ENABLE_RPATH)
include(AddInstallRPATHSupport)
add_install_rpath_support(BIN_DIRS "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}"
                          LIB_DIRS "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
                          DEPENDS ENABLE_RPATH
                          USE_LINK_PATH)

# Encourage user to specify a build type (e.g. Release, Debug, etc.), otherwise set it to Release.
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
      message(STATUS "Setting build type to 'Release' as none was specified.")
      set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE "Release")
  endif()
endif()

find_package(Threads)

set(LIBRARY_TARGET_NAME ${PROJECT_NAME})

add_library(${LIBRARY_TARGET_NAME} zupply.cpp ${CMAKE_CURRENT_SOURCE_DIR}/include/levi/external/zupply.h)
add_library(levi::${LIBRARY_TARGET_NAME} ALIAS ${LIBRARY_TARGET_NAME})
target_include_directories(${LIBRARY_TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories(${LIBRARY_TARGET_NAME} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(${LIBRARY_TARGET_NAME} PRIVATE Threads::Threads)

set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES PUBLIC_HEADER include/levi/external/zupply.h)

# install
install(TARGETS ${LIBRARY_TARGET_NAME}
        EXPORT levi
        COMPONENT runtime
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/levi/external)

if (NOT TARGET uninstall)
    include(AddUninstallTarget)
endif()
