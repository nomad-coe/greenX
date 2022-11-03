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
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_include_directories(LibGXAC PUBLIC src/)

target_sources(LibGXAC PRIVATE
        src/pade_approximant.f90
        )

target_link_libraries(LibGXAC GXCommon)

# -----------------------------------------------
# Library Installation
# -----------------------------------------------
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(TARGETS LibGXAC
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION include)

# -----------------------------------------------
# Unit Testing Set-up
# -----------------------------------------------
# Set name of test sub-directory in the build directory
set(TEST_TARGET_DIR "analytic-continuation")

# Create a directory in the build folder to place generated test drivers
set(UNIT_TEST_DIR ${PROJECT_BINARY_DIR}/unit_tests/${TEST_TARGET_DIR})
file(MAKE_DIRECTORY ${UNIT_TEST_DIR})
message("-- Analytic continuation unit tests written to: ${UNIT_TEST_DIR}")

# TODO(Alex)
# Everything from here can be placed into a function

# Generate test program from module
add_custom_command(
        OUTPUT ${UNIT_TEST_DIR}/test_pade_approximant_driver.f90
        COMMAND ${ZOFU_DRIVER} "${CMAKE_CURRENT_SOURCE_DIR}/src/test_pade_approximant.f90" "${UNIT_TEST_DIR}/test_pade_approximant_driver.f90"
        COMMENT "Generating ${UNIT_TEST_DIR}/test_pade_approximant_driver.f90"
    )

# Create the test driver executable
add_executable(test_pade_approximant)

# Set directory in which unit tests are built in
set_target_properties(test_pade_approximant
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${UNIT_TEST_DIR}"
        INCLUDE_DIRECTORIES "${ZOFU_INCLUDE_PATH}"
        )

target_sources(test_pade_approximant PRIVATE
        # Module
        src/test_pade_approximant.f90
        # Trivial, aut-generated program to execute module
        ${UNIT_TEST_DIR}/test_pade_approximant_driver.f90)

# Ensure our library gets compiled if one attempts to build the unit test executable
add_dependencies(test_pade_approximant LibGXAC)

target_link_libraries(test_pade_approximant ${LibZofu} LibGXAC)

add_test(NAME UNITTEST_pade_approximant
        COMMAND ${UNIT_TEST_DIR}/test_pade_approximant)
