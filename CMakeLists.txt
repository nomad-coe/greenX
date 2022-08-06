# GreenX Top-level CMake
#
# CMake Variable Descriptions
# PROJECT_BINARY_DIR:  Build directory of the most-recent project() command
# CMAKE_BINARY_DIR: Full path to the top level of the current CMake build tree
#
# A Buccheri.

cmake_minimum_required(VERSION 3.0.2)

project(greenX
    LANGUAGES Fortran
    VERSION 0.0.1
)

# Set folder structure for build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
# Note, this will put modules of each library in <BUILD_DIR>/modules
# We should consider if this is the desired behaviour
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
# Note, this will not find MKL or openBLAS
find_package(BLAS REQUIRED)

message('BLAS_FOUND: ${BLAS_FOUND}')
message('BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}')
message('BLAS_LIBRARIES: ${BLAS_LIBRARIES}')
message('BLAS95_FOUND: ${BLAS95_FOUND}')

find_package(LAPACK REQUIRED)
message('LAPACK_FOUND: ${LAPACK_FOUND}')
message('LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}')
message('LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}')

# Python required for application testing
include(cmake/python3.cmake)

# Enable ctest
enable_testing()

# Our library directories
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
add_subdirectory(GX-common)
add_subdirectory(GX-TimeFrequency)

# Set installation location
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(default_install_prefix "${PROJECT_SOURCE_DIR}/install")        
    set(CMAKE_INSTALL_PREFIX ${default_install_prefix}
            CACHE STRING "Choose the installation directory. Default location is ${default_install_prefix}"
            FORCE)
endif()
