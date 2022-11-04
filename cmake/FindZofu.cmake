# Find Zofu Unit Testing Framework
#
# Returns:
# ----------
#  LibZofu                Zofu library, statically compiled.
#  ZOFU_INCLUDE_PATH      Path to include folder.
#  ZOFU_DRIVER            Zofu driver binary, prepended by the full path to it.
#
# Notes on Finding Zofu
# ---------------------
# Option 1
#  Manually build in location of one's choosing, then pass ZOFU_PATH to CMake
# Option 2
#  If LibZofu is not found, pull, build and install.
#
# Zofu does not provide a package config file, such as zofu-config.cmake,
# therefore one has to the variables manually.
#
# Notes On Zofu Build and ExternalProject_Add
# -----------------------------------------------
# * On Mac, -DENABLE_MPI=ON did not work.
#   I had to add quotes to "${MPI_Fortran_INCLUDE_DIRS}" which appears in the
#   set_target_properties(...) calls, however the build still experienced problems.
#
# * ExternalProject_Add argument options should NOT be in quotes.
# * Use `add_library` to define the library in CMake. find_package won't necessarily
#   work because there's no guarantee that the library has already been build at
#   CMake configuration time.

# Pass path to Zofu install location to CMake
set(ZOFU_PATH "" CACHE STRING "Location of Zofu unit-testing library")

if (ZOFU_PATH)
    find_library(LibZofu NAME "libzofu" "zofu" HINTS "${ZOFU_PATH}/lib")

    # All targets get this added to their include path.
    # Note, the default is `finclude` so this will break if instructions are not followed
    set(ZOFU_INCLUDE_PATH ${ZOFU_PATH}/include)
    include_directories(${ZOFU_INCLUDE_PATH})

    # Program that generates a unit test driver given a test module.
    set(ZOFU_DRIVER ${ZOFU_PATH}/bin/zofu-driver)
endif()

# NOTE(ALEX)
# This all works, however line 72-73 in unit_test_functions.cmake will then require
# `target_link_libraries(${FUNC_TEST_NAME} LibZofu ${FUNC_REQUIRED_LIBS})`
# rather than
# `target_link_libraries(${FUNC_TEST_NAME} ${LibZofu} ${FUNC_REQUIRED_LIBS})`
# if ExternalProject_Add is used, and I do not immediately know how to reconcile this
# Hence, I've commented it out and provided build instructions for Zofu

# If the library is not found (either not installed, or the path is wrong) clone, build and install.
#if (NOT LibZofu)
#    if (NOT ZOFU_PATH STREQUAL "")
#        message("-- ZOFU_PATH was set but Zofu library was not found at that location.")
#    endif()
#
#    message("-- Cloning and installing Zofu")
#    set(ZOFU_PATH "${CMAKE_SOURCE_DIR}/external/zofu/install")
#
#    include(ExternalProject)
#
#    # Note, this will run at build time, NOT configure time
#    ExternalProject_Add(INTERAL_ZOFU
#            GIT_REPOSITORY https://github.com/acroucher/zofu     # Repo https
#            SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/zofu     # Location to clone to
#            GIT_SHALLOW    TRUE                                  # git clone --depth 1 to avoid downloading the whole history
#            GIT_PROGRESS   TRUE                                  # Report progress of git clone. More verbose CMake output
#            BUILD_ALWAYS   TRUE
#            CMAKE_ARGS -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=${ZOFU_PATH} -DZOFU_FORTRAN_MODULE_INSTALL_DIR:PATH=include
#            BUILD_COMMAND make
#            INSTALL_COMMAND make install
#            )
#
#    set(ZOFU_INCLUDE_PATH ${ZOFU_PATH}/include)
#    include_directories(${ZOFU_INCLUDE_PATH})
#
#    set(ZOFU_DRIVER ${ZOFU_PATH}/bin/zofu-driver)
#
#    add_library(LibZofu STATIC IMPORTED)
#    add_dependencies(LibZofu INTERAL_ZOFU)
#
#    set_target_properties(LibZofu PROPERTIES
#            IMPORTED_LOCATION "${ZOFU_PATH}/lib/libzofu.a"
#            INTERFACE_INCLUDE_DIRECTORIES "${ZOFU_INCLUDE_PATH}"
#            )
#endif ()


if (LibZofu)
    message("-- Found LibZofu ${LibZofu}")
    message("-- LibZofu's module path: ${ZOFU_INCLUDE_PATH}")
else()
    message("-- LibZofu not built at configure time")
endif()
