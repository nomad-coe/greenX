# Define an application test

function(add_app_test)
    # Create an application test comprised of a fortran executable and a
    # python driver (to provide input, run the exe and assert the result).
    #
    # Expectations
    #  * Naming convention follows: TEST_NAME.f90 TEST_NAME.py
    #    defining the app test test_TEST_NAME.
    #  * Tests are defined in the test/ folder, relative to whichever library
    #   they are written for.
    #
    # https://cmake.org/cmake/help/latest/command/cmake_parse_arguments.html

    # Define the CMake function signature keywords and their types
    # Of the general form: set(type KEYWORD)
    set(oneValueArgs TEST_NAME TEST_TARGET_DIR)     # Single-value options
    set(multiValueArgs LIBS_FOR_TESTING)            # Multi-value options: Multiple arguments or list/s
    cmake_parse_arguments(FUNC                      # Prefix for all function arguments within function body
            "${options}"                            # Assign the binary options for the function
            "${oneValueArgs}"                       # Assign the single-value options for the function
            "${multiValueArgs}"                     # Assign the multi-value options for the function
            ${ARGV})                                # ${ARGN} or ${ARGV}. (I think) ${ARGV} means accept a variable
                                                    # number of arguments, which one want for a list of no fixed size

    # Set target name
    add_executable(${TEST_NAME})

    # Set binary name of target
    set_target_properties(${TEST_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_NAME "${TEST_NAME}.exe")

    # Define source that comprise the binary
    target_sources(${TEST_NAME}
            PRIVATE
            test/${TEST_NAME}.f90
            )

    # Libraries that the binary links to
    target_link_libraries(${TEST_NAME}
            PUBLIC
            ${LIBS_FOR_TESTING}
            )

    # Build location of the binary
    set_target_properties(${TEST_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})

    # Copy .py test to the `<build>/test/${WORKING_TEST_DIR}` directory, such that one can
    #  `<build>/test/${WORKING_TEST_DIR} && pytest -s`
    # or
    #  `<build>/test && pytest -s`, or
    add_custom_command(
            TARGET LibGXMiniMax POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            # Test source relative to the time-frequency folder
            # where CMAKE_CURRENT_SOURCE_DIR => CMakeLists.txt on this level
            ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_NAME}.py
            # Location to copy the test to
            ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}/${TEST_NAME}.py)

    # Add test to ctest
    add_test(
            NAME ${TEST_NAME}
            COMMAND pytest -s
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/${TEST_TARGET_DIR}
    )

endfunction()
