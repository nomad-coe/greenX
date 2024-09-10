# ***************************************************************************************************
#  Copyright (C) 2020-2024 GreenX library
#  This file is distributed under the terms of the APACHE2 License.
#
# ***************************************************************************************************

# -----------------------------------------------
# CMakeLists for Analytic Continuation Library
# -----------------------------------------------
add_library(LibGXAC "")

set_target_properties(LibGXAC
  PROPERTIES
  VERSION 0.0.1
  SOVERSION 0.0.1
  LIBRARY_OUTPUT_NAME gx_ac
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY})

target_include_directories(LibGXAC PUBLIC src/)
if(GMPXX_FOUND)
  add_definitions(-DGMPXX_FOUND)
  target_sources(LibGXAC PRIVATE src/pade_approximant.f90 src/ComplexGMP.cpp src/Symmetry_pade.cpp src/pade_mp.cpp api/gx_ac.F90)
  target_include_directories(LibGXAC PRIVATE ${GMPXX_INCLUDE_DIRS})
  target_link_libraries(LibGXAC GXCommon ${GMPXX_LIBRARIES} gmp)
else()
  target_sources(LibGXAC PRIVATE src/pade_approximant.f90 api/gx_ac.F90)
  target_link_libraries(LibGXAC GXCommon)
endif()

# -----------------------------------------------
# Library Installation
# -----------------------------------------------
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(TARGETS LibGXAC ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}  DESTINATION include)



# -----------------------------------------------
# Application Testing Set-Up
# -----------------------------------------------
# Include cmake custom function
include(../cmake/testFunctions.cmake)

# Set name of test sub-directory in the build directory
set(TEST_TARGET_DIR "analytic_continuation")

target_include_directories(LibGXAC PUBLIC test/)

# Set pytest conftest for Localized basis library tests
add_custom_command(
	TARGET LibGXAC POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/test/conftest.py
        ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/conftest.py)


# -----------------------------------------------
# Application Tests
# -----------------------------------------------

# Add a test target
add_executable(test_gx_analytic_continuation)

# Set binary name
set_target_properties(test_gx_analytic_continuation
        PROPERTIES
        RUNTIME_OUTPUT_NAME test_gx_analytic_continuation.exe)

# Define source that comprise the binary
target_sources(test_gx_analytic_continuation
        PRIVATE
        test/test_gx_analytic_continuation.f90
        )

# Libraries that the binary links to
target_link_libraries(test_gx_analytic_continuation
        PUBLIC
	    LibGXAC
        )

# Build location of the binary
# TODO(Alex) Consider move this to test/, to sit with the python drivers
set_target_properties(test_gx_analytic_continuation
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})

# Copy .py test to the `build/test` directory, such that one can run pytest there
# where CMAKE_CURRENT_SOURCE_DIR => CMakeLists.txt on this level
add_custom_command(
	TARGET LibGXAC POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency (this) folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_analytic_continuation_no-greedy_64bit.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_analytic_continuation_no-greedy_64bit.py)
add_custom_command(
	TARGET LibGXAC POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency (this) folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_analytic_continuation_greedy_64bit.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_analytic_continuation_greedy_64bit.py)
add_custom_command(
	TARGET LibGXAC POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency (this) folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_analytic_continuation_no-greedy_128bit.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_analytic_continuation_no-greedy_128bit.py)
add_custom_command(
	TARGET LibGXAC POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency (this) folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_analytic_continuation_greedy_128bit.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/test_gx_analytic_continuation_greedy_128bit.py)

# Add test to ctest
add_test(
        NAME test_gx_analytic_continuation_no-greedy_64bit
        COMMAND pytest -s test_gx_analytic_continuation_no-greedy_64bit.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)
add_test(
        NAME test_gx_analytic_continuation_greedy_64bit
        COMMAND pytest -s test_gx_analytic_continuation_greedy_64bit.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)
add_test(
        NAME test_gx_analytic_continuation_no-greedy_128bit
        COMMAND pytest -s test_gx_analytic_continuation_no-greedy_128bit.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)
add_test(
        NAME test_gx_analytic_continuation_greedy_128bit
        COMMAND pytest -s test_gx_analytic_continuation_greedy_128bit.py --build-dir ${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
)