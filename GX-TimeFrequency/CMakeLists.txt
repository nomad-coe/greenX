# **************************************************************************************************
#  Copyright (C) 2020-2023 GreenX library
#  This file is distributed under the terms of the APACHE2 License.
#
# **************************************************************************************************

# CMakeLists for TimeFrequency Library (MiniMax)

# GreenXMiniMax Library
add_library(LibGXMiniMax "")

# Set properties of the library
set_target_properties(LibGXMiniMax
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME gx_minimax
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Add include directories to the LibGXMiniMax
#target_include_directories(LibGXMiniMax PUBLIC src/ api/ utilities/)
target_include_directories(LibGXMiniMax
    PUBLIC
        # Build-time includes (point to actual source dirs)
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utilities>

        # Install-time includes (point to the location where headers/modules land when installed)
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_Fortran_MODULES}/TimeFrequency/modules>
)
# Define source that comprise LibGXMiniMax
target_sources(LibGXMiniMax PRIVATE
        src/gx_common.h
	src/minimax_utils.F90
        src/minimax_tau.F90
        src/minimax_omega.F90
        src/minimax_grids.F90
        api/api_utilities.f90
        api/gx_minimax.f90
        )

# Link external libraries and lower-level GX lib dependency
target_link_libraries(LibGXMiniMax GXCommon ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})


# -----------------------------------------------
# Build utility program for grid tabulation
# -----------------------------------------------
add_executable(GXTabulateGrids)

set_target_properties(GXTabulateGrids
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY}
        RUNTIME_OUTPUT_NAME gx_tabulate_grids.exe)

target_sources(GXTabulateGrids
        PRIVATE
        utilities/gx_tabulate_minimax.F90
        src/gx_common.h
        )

target_link_libraries(GXTabulateGrids
        PRIVATE
        GXCommon
        LibGXMiniMax
        )

# -----------------------------------------------
# Library Installation
# -----------------------------------------------
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(TARGETS LibGXMiniMax
        EXPORT greenXTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION ${CMAKE_INSTALL_Fortran_MODULES}/TimeFrequency)

# Install `Build utility` program
install(TARGETS GXTabulateGrids
        DESTINATION bin)

# -----------------------------------------------
# Application Testing Set-Up
# -----------------------------------------------
# Include cmake custom function
include(../cmake/testFunctions.cmake)

# Set name of test sub-directory in the build directory
set(TEST_TARGET_DIR "time-frequency")

target_include_directories(LibGXMiniMax 
         PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test>
)

# Set pytest conftest for Time-Frequency library tests
add_custom_command(
        TARGET LibGXMiniMax POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/conftest.py
        ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/conftest.py)

# -----------------------------------------------
# Application Tests
# -----------------------------------------------

# Add a test target
add_executable(test_gx_minimax_grid)

# Set binary name
set_target_properties(test_gx_minimax_grid
        PROPERTIES
        RUNTIME_OUTPUT_NAME test_gx_minimax_grid.exe)

# Define source that comprise the binary
target_sources(test_gx_minimax_grid
        PRIVATE
        test/test_gx_minimax_grid.f90
        )

# Libraries that the binary links to
target_link_libraries(test_gx_minimax_grid
        PUBLIC
        LibGXMiniMax
        )

# Build location of the binary
# TODO(Alex) Consider move this to test/, to sit with the python drivers
set_target_properties(test_gx_minimax_grid
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})

# Copy .py test to the `build/test` directory, such that one can run pytest there
# where CMAKE_CURRENT_SOURCE_DIR => CMakeLists.txt on this level
add_custom_command(
        TARGET LibGXMiniMax POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency (this) folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_minimax_grid.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_minimax_grid.py)

# Note.
# The command below would build the test in <BUILD_DIR>/<LIB_NAME>/test:
#                #  Location to copy the test to
#                ${CMAKE_CURRENT_BINARY_DIR}/test/test_gx_minimax_grid.py)
# whereas I am currently opting to put tests in <BUILD_DIR>/tests/<LIB_NAME>

# Add test to ctest
add_test(
        NAME test_gx_minimax_grid
        COMMAND pytest -s test_gx_minimax_grid.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)

# Add a second test, using the `tabulate utility` binary
add_custom_command(
        TARGET LibGXMiniMax POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_tabulate_minimax.py
        ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_tabulate_minimax.py)

add_test(
        NAME test_gx_tabulate_minimax
        COMMAND pytest -s test_gx_tabulate_minimax.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)


# TODO For developers.
# To make the addition of further tests more compact, one can use this function call
# For example, for the test above:
#
#```
#   set(LIBS_FOR_TESTING)
#   list(APPEND LIBS_FOR_TESTING LibGXMiniMax)
#   set(TEST_NAME "test_gx_minimax_grid")
#   add_app_test(TEST_NAME TEST_TARGET_DIR LIBS_FOR_TESTING)
#```
