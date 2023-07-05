# Python helper functions

function(find_python_module module)
    # Find a python package.
    # Returns ${module}_FOUND
    #
    # Examples:
    # Find an optional python package, with any version number:
    #  `find_python_module(pygreenx)`
    #
    # Find a required python package, with a specific version number:
    #   `find_python_module(pygreenx REQUIRED VERSION 1.1.1)`
    #
    # Note, https://cmake.org/cmake/help/latest/module/FindPython.html
    # gives variables as Python_FOUND, etc, however one actually may need to specify
    # Python3_FOUND, etc.

    # Function arguments
    set(options REQUIRED)
    set(oneValueArgs VERSION)

    # Parse function arguments, prepended with MY_
    cmake_parse_arguments(MY "${options}" "${oneValueArgs}" "" ${ARGN})

    if(NOT DEFINED Python3_EXECUTABLE)
        message(FATAL_ERROR "No python interpreter has been found by CMake.")
    endif()

    # Shell command to query python package availability
    execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import ${module}; print(${module}.__version__)"
            RESULT_VARIABLE CHECK_RESULT
            OUTPUT_VARIABLE MODULE_VERSION
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Annoyingly, one cannot assign error type for passing to message, so the code is duplicated
    # Errors
    if(MY_REQUIRED)
        # Restrictions on version
        if(DEFINED MY_VERSION)
            if(CHECK_RESULT EQUAL 0 AND MODULE_VERSION VERSION_GREATER_EQUAL MY_VERSION)
                message(STATUS "Found ${module}: ${MODULE_VERSION}")
                set(${module}_FOUND TRUE PARENT_SCOPE)
            elseif(CHECK_RESULT EQUAL 0 AND NOT MODULE_VERSION VERSION_GREATER_EQUAL MY_VERSION)
                message(FATAL_ERROR "${module} version is not sufficient. Require ${MY_VERSION} but found ${MODULE_VERSION}")
            else()
                message(FATAL_ERROR "${module} was not found.")
            endif()

        # No restrictions on version
        else()
            if(CHECK_RESULT EQUAL 0)
                message(STATUS "Found ${module}: ${MODULE_VERSION}")
                set(${module}_FOUND TRUE PARENT_SCOPE)
            else()
                message(FATAL_ERROR "${module} was not found.")
            endif()
        endif()

    # Warnings
    else()
        # Restrictions on version
        if(DEFINED MY_VERSION)
            if(CHECK_RESULT EQUAL 0 AND MODULE_VERSION VERSION_GREATER_EQUAL MY_VERSION)
                message(STATUS "Found ${module}: ${MODULE_VERSION}")
                set(${module}_FOUND TRUE PARENT_SCOPE)
            elseif(CHECK_RESULT EQUAL 0 AND NOT MODULE_VERSION VERSION_GREATER_EQUAL MY_VERSION)
                message(WARNING "${module} version is not sufficient. Require ${MY_VERSION} but found ${MODULE_VERSION}")
                set(${module}_FOUND FALSE PARENT_SCOPE)
            else()
                message(WARNING "${module} was not found.")
                set(${module}_FOUND FALSE PARENT_SCOPE)
            endif()

        # No restrictions on version
        else()
            # Module is installed. No restrictions on version
            if(CHECK_RESULT EQUAL 0)
                message(STATUS "Found ${module}: ${MODULE_VERSION}")
                set(${module}_FOUND TRUE PARENT_SCOPE)
            else()
                message(WARNING "${module} was not found.")
                set(${module}_FOUND FALSE PARENT_SCOPE)
            endif()
        endif()
    endif()

endfunction()


function(find_pythonhome)
    if ($ENV{PYTHONHOME})
        set(_PYTHONHOME $ENV{PYTHONHOME})
    else()
        execute_process(
                COMMAND "${Python3_EXECUTABLE}" "-c"
                "import sys; print(sys.prefix + ':' + sys.exec_prefix)"
                RESULT_VARIABLE _PYTHONHOME_FAILED
                OUTPUT_VARIABLE _PYTHONHOME
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(_PYTHONHOME_FAILED)
            message(FATAL_ERROR "Could not determine PYTHONHOME. Error:" ${_PYTHONHOME_FAILED})
        endif()
    endif()
    message("    Found PYTHONHOME: " ${_PYTHONHOME})
    set(PYTHONHOME ${_PYTHONHOME} PARENT_SCOPE)
endfunction(find_pythonhome)
