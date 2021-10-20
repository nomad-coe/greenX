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

# External libraries
# Python required for application testing
find_package(Python3 3.7 COMPONENTS Interpreter Development)
if(Python3_FOUND)
    message("-- Python 3 interpreter version: " ${Python3_VERSION})
else()
    message("-- Python 3 interpreter not found")
endif()

# Enable ctest
enable_testing()

# Our libraries 
add_subdirectory(GX-TimeFrequency)
