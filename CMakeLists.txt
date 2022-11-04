# This file is distributed under the terms of the APACHE2 License.
# GreenX Top-level CMakeLists.txt
#
# CMake Variable Descriptions
# PROJECT_BINARY_DIR:  Build directory of the most-recent project() command
# CMAKE_BINARY_DIR: Full path to the top level of the current CMake build tree

cmake_minimum_required(VERSION 3.15.0)

project(greenX
    LANGUAGES Fortran
    VERSION 0.0.1
)

# Set folder structure for build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
# Note, this will put modules of each library in <BUILD_DIR>/modules
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Define GNU standard installation directories
include(GNUInstallDirs)

# CMake module directory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Defines build type
include(cmake/StandardProjectSettings.cmake)

# Compiler flags
include(cmake/CompilerFlags.cmake)

# Compiler warnings
include(cmake/CompilerWarnings.cmake)

# External libraries
# Note, for Blas and Lapack, this is also able to find vendor-specific versions:
# https://cmake.org/cmake/help/v3.15/module/FindBLAS.html
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Python required for application testing
include(cmake/python3.cmake)
find_python_module(numpy)
find_python_module(pytest)

# Enable ctest
enable_testing()

# Optional unit testing lib
option(ENABLE_UNITTESTS "Enable GreenX Unit Testing" OFF)
if (${ENABLE_UNITTESTS})
    # Function for generating unit test executables (assumes Zofu)
    include(cmake/unit_test_functions.cmake)
    # Calling find_package is equivalent to: include(cmake/FindZofu.cmake)
    find_package(zofu REQUIRED)
endif()

# Our library directories
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
add_subdirectory(GX-AnalyticContinuation)
add_subdirectory(GX-common)
add_subdirectory(GX-TimeFrequency)

# Documentation
option(ENABLE_DOCS "Enable documentation" ON)
if (${ENABLE_DOCS})
    find_program(DOXYGEN doxygen REQUIRED)
    message("-- Doxygen documentation support enabled.")
    set(doxy_output "documentation")
    file(MAKE_DIRECTORY ${doxy_output})
    add_custom_target(docs COMMAND ${DOXYGEN} ${PROJECT_SOURCE_DIR}/Doxyfile)
endif ()

# Set installation location
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(default_install_prefix "${PROJECT_SOURCE_DIR}/install")        
    set(CMAKE_INSTALL_PREFIX ${default_install_prefix}
            CACHE STRING "Choose the installation directory. Default location is ${default_install_prefix}"
            FORCE)
endif()
