# Find Zofu Unit Testing Framework
#
# Returns:
# ----------
#  LibZofu
#  ?
#  ?
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


# First approach. Pass path to Zofu install location to CMake
set(ZOFU_PATH "" CACHE STRING "Location of Zofu unit-testing library")
find_library(LibZofu NAME "libzofu" "zofu" HINTS "${ZOFU_PATH}/lib")

# More sensible include name, suggested in the docs and used when caling ExternalProject_Add
# include_directories(${ZOFU_PATH}/finclude)
set(ZOFU_INCLUDE_PATH ${ZOFU_PATH}/include)
include_directories(${ZOFU_INCLUDE_PATH})

# Program that generates a unit test driver given a test module.
set(ZOFU_DRIVER ${ZOFU_PATH}/bin/zofu-driver)

# If the library is not found (either not installed, or the path is wrong)
# clone, build and install.
if (NOT LibZofu)
    if (NOT ZOFU_PATH STREQUAL "")
        message("-- ZOFU_PATH was set but Zofu library was not found at that location.")
    endif()
    message("-- Cloning and installing Zofu")
    set(ZOFU_PATH "${CMAKE_SOURCE_DIR}/external/zofu/install")

    include(ExternalProject)


#    ExternalProject_Add(INTERAL_ZOFU
#            GIT_REPOSITORY https://github.com/acroucher/zofu     # Repo https
#            GIT_SHALLOW    TRUE                                  # git clone --depth 1 to avoid downloading the whole history
#            GIT_PROGRESS   TRUE                                  # Report progress of git clone. More verbose CMake output
#            SOURCE_DIR     "${CMAKE_SOURCE_DIR}/external/zofu/"  # Location to clone to
#            BUILD_ALWAYS   FALSE
#            BUILD_IN_SOURCE 1
#            CMAKE_ARGS -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=${ZOFU_PATH} -DZOFU_FORTRAN_MODULE_INSTALL_DIR:PATH=include
#            BUILD_COMMAND make
#            INSTALL_COMMAND make install
#            )

    add_library(LibZofu STATIC IMPORTED)
    set_target_properties(LibZofu PROPERTIES IMPORTED_LOCATION ${ZOFU_PATH}/lib/libzofu.a)
    # All targets get this added to their include path
    include_directories(${ZOFU_INCLUDE_PATH})
endif ()

if (LibZofu)
    message("-- Found LibZofu ${LibZofu}")
    message("-- LibZofu's module path: ${ZOFU_INCLUDE_PATH}")

else()
    message("-- LibZofu not found")
    # TODO(Alex) Get this to return an error
endif()
