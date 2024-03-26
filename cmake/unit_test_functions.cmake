# Create unit test binary from a module that uses the Zofu unit testing framework
# LibZofu and ZOFU_INCLUDE_PATH are implicitly available as globals.
#
# Conventions
# --------------
# The test module should be named ${TEST_NAME}.f90
# The test binary will have the name ${TEST_NAME}
# A unit test with be created with the CMake variable name UNITTEST_${TEST_NAME}
# ${TARGET_TEST_DIR} should be somewhere senible in the build directory, for
# example build/unit_tests/library_name
# The test driver - the auto-generated program that calls the module ${TEST_NAME}.f90 -
# is placed in the same location as the test binary: ${TARGET_TEST_DIR}
#
# Arguments
# --------------
# TARGET_TEST_DIR: Build folder subdirectory in which to save the driver.f90
#   and corresponding unit test binary.
# TEST_NAME: Name of the test module (excluding the .f90 extension)
#   This will also be used to name the unit test
# REQUIRED_LIBS: List of libaries (excluding Zofu) that the unit tests depend on.
#
# Notes
# ------------
# One must link to LibZofu rather than ${LibZofu} when externalProject_add used.
# Hence externalProject_add is commented-out in Findzofu.cmake.
# One should probably look into the nuance of linking to LIB vs ${LIB}
# https://cmake.org/cmake/help/latest/guide/importing-exporting/index.html#importing-libraries

function(create_unit_test_executable)
    # Function arguments:
    set(oneValueArgs TARGET_TEST_DIR TEST_NAME)     # Single-value options
    set(multiValueArgs REQUIRED_LIBS)               # Multi-value options: Multiple arguments or list/s

    # Parse function arguments and prepend prefix
    # This should be present in all functions, with the same signature.
    # The prefix gives variables a local scope.
    cmake_parse_arguments(FUNC                      # Prefix for all function arguments within function body
            "${options}"                            # Assign the binary options for the function
            "${oneValueArgs}"                       # Assign the single-value options for the function
            "${multiValueArgs}"                     # Assign the multi-value options for the function
             ${ARGV})                               # ${ARGN} or ${ARGV}. (I think) ${ARGV} means accept a variable

    # Generate test program .f90 from module
    add_custom_command(
            OUTPUT "${FUNC_TARGET_TEST_DIR}/${FUNC_TEST_NAME}_driver.f90"
            # For example, ./zofu-driver path/to/module.f90 target/location/module_driver.f90
            COMMAND ${ZOFU_DRIVER} "${CMAKE_CURRENT_SOURCE_DIR}/test/${FUNC_TEST_NAME}.f90" "${FUNC_TARGET_TEST_DIR}/${FUNC_TEST_NAME}_driver.f90"
            COMMENT "Generating ${FUNC_TARGET_TEST_DIR}/${FUNC_TEST_NAME}_driver.f90"
    )

    # Create the test executable
    add_executable(${FUNC_TEST_NAME})

    # Set directory in which unit tests are built in
    set_target_properties(${FUNC_TEST_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${FUNC_TARGET_TEST_DIR}"
            INCLUDE_DIRECTORIES "${ZOFU_INCLUDE_PATH}"
            )

    # Specify source code that the target will depend on
    target_sources(${FUNC_TEST_NAME} PRIVATE
            # Test module
            test/${FUNC_TEST_NAME}.f90
            # Trivial, aut-generated program to execute module
            "${FUNC_TARGET_TEST_DIR}/${FUNC_TEST_NAME}_driver.f90")

    # Ensure the library we're testing gets compiled if one attempts to build
    # the corresponding unit test executable
    add_dependencies(${FUNC_TEST_NAME} ${FUNC_REQUIRED_LIBS})

    # Link libraries to the test
    target_link_libraries(${FUNC_TEST_NAME} ${LibZofu} ${FUNC_REQUIRED_LIBS})

    # Add test to ctest
    add_test(NAME UNITTEST_${FUNC_TEST_NAME}
             COMMAND ${FUNC_TARGET_TEST_DIR}/${FUNC_TEST_NAME})

endfunction()
