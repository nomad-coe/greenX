# GreenX Top-level CMake 
#
# CMake Variable Descriptions 
# PROJECT_BINARY_DIR:  Build directory of the most-recent project() command
# CMAKE_BINARY_DIR: Full path to the top level of the current CMake build tree

cmake_minimum_required(VERSION 3.0.2)

project(greenX 
    LANGUAGES Fortran
    VERSION 0.0.1
)

# Set folder structure for build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
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

# Python required for application testing
include(cmake/python3.cmake)

# Enable ctest
enable_testing()

# Our libraries
add_subdirectory(GX-common)
add_subdirectory(GX-TimeFrequency)
