# **************************************************************************************************
#  Copyright (C) 2020-2023 GreenX library
#  This file is distributed under the terms of the APACHE2 License.
#
# **************************************************************************************************

# CMakeLists for Localized basis Library (NAOs)

# GreenXLBasis Library

# CMakeLists for Localized Basis component of the GreenX Library

# GreenX Localized Basis Library
add_library(LibGXLBasis "")

# Set properties of the library
set_target_properties(LibGXLBasis
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME gx_localized_basis
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Add include directories to the LibGXMiniMax
target_include_directories(LibGXLBasis PUBLIC src/ api/)

# Define source that comprise LibGXLBasis
target_sources(LibGXLBasis PRIVATE
	src/localized_basis.f90
        src/localized_basis_types.f90
	src/localized_basis_environments.f90
	api/gx_loc_basis.f90
        )

# Link external libraries and lower-level GX lib dependency
target_link_libraries(LibGXLBasis GXCommon ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# -----------------------------------------------
# Library Installation
# -----------------------------------------------
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(TARGETS LibGXLBasis
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION include)

# -----------------------------------------------
# Application Testing Set-Up
# -----------------------------------------------
# Include cmake custom function
include(../cmake/testFunctions.cmake)

# Set name of test sub-directory in the build directory
set(TEST_TARGET_DIR "localized_basis")

target_include_directories(LibGXLBasis PUBLIC test/)

# Set pytest conftest for Localized basis library tests
add_custom_command(
	TARGET LibGXLBasis POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/conftest.py
        ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/conftest.py)

# Set input data (density) for Localized basis library tests
add_custom_command(
        TARGET LibGXLBasis POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/density.dat
        ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/density.dat)

# Set input data (density) for Localized basis library tests
add_custom_command(
        TARGET LibGXLBasis POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/overlap.dat
        ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/overlap.dat)


# -----------------------------------------------
# Application Tests
# -----------------------------------------------

# Add a test target
add_executable(test_gx_localized_basis)

# Set binary name
set_target_properties(test_gx_localized_basis
        PROPERTIES
        RUNTIME_OUTPUT_NAME test_gx_localized_basis.exe)

# Define source that comprise the binary
target_sources(test_gx_localized_basis
        PRIVATE
        test/test_gx_localized_basis.f90
        )

# Libraries that the binary links to
target_link_libraries(test_gx_localized_basis
        PUBLIC
	LibGXLBasis
        )

# Build location of the binary
# TODO(Alex) Consider move this to test/, to sit with the python drivers
set_target_properties(test_gx_localized_basis
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})

# Copy .py test to the `build/test` directory, such that one can run pytest there
# where CMAKE_CURRENT_SOURCE_DIR => CMakeLists.txt on this level
add_custom_command(
	TARGET LibGXLBasis POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency (this) folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_localized_basis.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_localized_basis.py)

# Note.
# The command below would build the test in <BUILD_DIR>/<LIB_NAME>/test:
#                #  Location to copy the test to
#                ${CMAKE_CURRENT_BINARY_DIR}/test/test_gx_localized_basis.py)
# whereas I am currently opting to put tests in <BUILD_DIR>/tests/<LIB_NAME>

# Add test to ctest
add_test(
        NAME test_gx_localized_basis
        COMMAND pytest -s test_gx_localized_basis.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)

