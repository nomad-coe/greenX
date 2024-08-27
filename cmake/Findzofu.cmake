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

# if ZOFU_PATH given, try to find it
if (ZOFU_PATH)
  find_library(LibZofu NAME "libzofu" "zofu" HINTS "${ZOFU_PATH}/lib")
endif()

# if found, use it; else: build on the fly
if (LibZofu)
  message("-- Found LibZofu ${LibZofu}")

  # Check for multiple possible include directories
  if (EXISTS "${ZOFU_PATH}/include")
    set(ZOFU_INCLUDE_PATH "${ZOFU_PATH}/include")
  elseif (EXISTS "${ZOFU_PATH}/finstall/zofu")
    set(ZOFU_INCLUDE_PATH "${ZOFU_PATH}/finstall/zofu")
  else()
    message(FATAL_ERROR "No suitable Zofu include directory found.")
  endif()
  include_directories(${ZOFU_INCLUDE_PATH})

  # Program that generates a unit test driver given a test module.
  set(ZOFU_DRIVER ${ZOFU_PATH}/bin/zofu-driver)

  message("-- LibZofu's module path: ${ZOFU_INCLUDE_PATH}")
else()
  message("-- LibZofu not found")
  # Build zofu at compile time
  set(ZOFU_BUILD_ON_FLY TRUE) 
  set(ZOFU_PATH_ROOT "${CMAKE_BINARY_DIR}/external/zofu")

  ExternalProject_Add(
      zofu
      PREFIX ${ZOFU_PATH_ROOT}  # Directory to download and build the project
      GIT_REPOSITORY https://github.com/acroucher/zofu.git  # GitHub repository URL
      GIT_TAG master  # Specific branch, tag, or commit
      CMAKE_ARGS
          -DCMAKE_BUILD_TYPE=release
          -DCMAKE_INSTALL_PREFIX=${ZOFU_PATH_ROOT}/install
          -DZOFU_FORTRAN_MODULE_INSTALL_DIR:PATH=include
      UPDATE_COMMAND ""  # Leave this empty to avoid re-downloading on every build
      TEST_COMMAND ""  # No test step
  )

  # set the zofu library 
  set(LibZofu "${ZOFU_PATH_ROOT}/install/lib/libzofu.a")

  # set the zofu library include directory
  set(ZOFU_INCLUDE_PATH "${ZOFU_PATH_ROOT}/install/include")
  include_directories(${ZOFU_INCLUDE_PATH})

  # Program that generates a unit test driver given a test module.
  set(ZOFU_DRIVER ${ZOFU_PATH_ROOT}/install/bin/zofu-driver)

  message("-- LibZofu will be built in ${ZOFU_PATH_ROOT}")
  message("-- LibZofu's module path: ${ZOFU_INCLUDE_PATH}")
endif()
