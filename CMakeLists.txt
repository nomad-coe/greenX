# This file is distributed under the terms of the APACHE2 License.
# GreenX Top-level CMakeLists.txt
#
# CMake Variable Descriptions
# PROJECT_BINARY_DIR:  Build directory of the most-recent project() command
# CMAKE_BINARY_DIR: Full path to the top level of the current CMake build tree

cmake_minimum_required(VERSION 3.15.0)

project(greenX
  LANGUAGES Fortran
  VERSION 0.0.1
)

# Need for Intel MKL and for GNU GMP
enable_language(CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set folder structure for build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Define GNU standard installation directories
include(GNUInstallDirs)

# CMake module directory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Defines build type
include(cmake/StandardProjectSettings.cmake)

# Compiler flags
include(cmake/CompilerFlags.cmake)

# Compiler warnings
include(cmake/CompilerWarnings.cmake)

# Allows CMake to compile external projects and submodules.
include(ExternalProject)

# External libraries
# Note, for Blas and Lapack, this is also able to find vendor-specific versions:
# https://cmake.org/cmake/help/v3.15/module/FindBLAS.html
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

option(ENABLE_GREENX_CTEST "Enable GreenX CTest" ON)
if (${ENABLE_GREENX_CTEST})
   # Python required for application testing
   include(cmake/python3.cmake)
   find_package(Python3 3.7 COMPONENTS Interpreter Development)
   if (Python3_FOUND)
      message("-- Python 3 interpreter version: " ${Python3_VERSION})
   else()
      message("-- Python 3 interpreter not found")
   endif()

   # Enable ctest
   enable_testing()

   # Application test python module
   # Note, pygreenx dependencies should be handled pygreenx
   # pygreenx version specified in its setup.cfg
   find_python_module(pygreenx VERSION 0.0.1)
   if (NOT pygreenx_FOUND)
      message("-- pygreenx is required for application testing. To install it, cd <PROJECTROOT>/python and run `pip install -e .`")
   endif ()
endif()   

# Optional unit testing lib
option(ENABLE_GREENX_UNIT_TESTS "Enable GreenX Unit Testing" OFF)

# Our library directories
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
add_subdirectory(GX-common)
option(AC_COMPONENT "Enable the Analytical Continuation component" ON)
if (${AC_COMPONENT})
   # Optional multiple precision arithmetic
   option(ENABLE_GNU_GMP "Enable the GNU GMP library for multiple precision arithmetic" ON)
   if (${ENABLE_GNU_GMP})
      find_package(GMPXX)
   endif()
   add_subdirectory(GX-AnalyticContinuation)
endif()
option(MINIMAX_COMPONENT "Enable the minimax time-frequency grids component" ON)
if (${MINIMAX_COMPONENT})
   add_subdirectory(GX-TimeFrequency)
endif ()

option (LBASIS_COMPONENT "Enable the localized basis component" ON) 
if (${LBASIS_COMPONENT})
   add_subdirectory(GX-LocalizedBasis)
endif()

# Elements to compile the libraries inside submodules
option(COMPILE_SUBMODULES "Compile GreenX component contained inside submodules" OFF)
if (${COMPILE_SUBMODULES})
  # zofu is needed for IDieL regression tests
  if (${ENABLE_GREENX_UNIT_TESTS})
    # Function for generating unit test executables (assumes Zofu)
    include(cmake/unit_test_functions.cmake)
    # find ZOFU in $ZOFU_PATH or build on the fly
    include(cmake/Findzofu.cmake)
  endif()
  # Build IDieL
  set(IDieL_SOURCE "${PROJECT_SOURCE_DIR}/GX-q=0/IDieL/")
  set(IDieL_CMAKE_OPTIONS "" CACHE STRING "Other IDieL options to build IDieL")
  # Check if IDieL source directory exists
  if (EXISTS ${IDieL_SOURCE})
    ExternalProject_Add(IDieL
        SOURCE_DIR ${IDieL_SOURCE}
        PREFIX ${IDieL_SOURCE}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                   -DCMAKE_Fortran_COMPILER=gfortran
                   ${IDieL_CMAKE_OPTIONS}
        BUILD_COMMAND make -j4
        INSTALL_COMMAND make install
        TEST_COMMAND cd ${IDieL_SOURCE}/src/IDieL-build/test && ctest --output-on-failure -C $<CONFIGURATION>
    )
  # Now properly link our IDieL installation from the IDieL install directory to the common one from GreenX
  set(IDieL_INSTALL_DIR "${IDieL_SOURCE}/IDieL")
  install(FILES ${IDieL_INSTALL_DIR}/lib/libIDieL.so DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  install(DIRECTORY ${IDieL_INSTALL_DIR}/include/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/modules/
        FILES_MATCHING PATTERN "*d")
  else()
    message(FATAL_ERROR "IDieL source directory not found: ${IDieL_SOURCE}. You can obtain it by executing 'git submodule update --init'.")
  endif()

endif()

# Install libpaw
option(PAW_COMPONENT "Enable the PAW component" OFF)
if (${PAW_COMPONENT})
   set(LIBPAW_SOURCE "${PROJECT_SOURCE_DIR}/GX-PAW/")
   ExternalProject_Add(LIBPAW
      SOURCE_DIR ${LIBPAW_SOURCE}
      PREFIX ${LIBPAW_SOURCE}
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                 -DCMAKE_Fortran_COMPILER=gfortran
      BUILD_COMMAND make -j 1
      INSTALL_COMMAND make install)
   # Now properly link our LIBPAW installation from the LIBPAW install directory to the common one from GreenX
   set(LIBPAW_INSTALL_DIR "${LIBPAW_SOURCE}/libpaw-greenx")
   install(FILES ${LIBPAW_INSTALL_DIR}/lib/libabinit_common.a DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
   install(FILES ${LIBPAW_INSTALL_DIR}/lib/libpaw.so DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
   install(DIRECTORY ${LIBPAW_INSTALL_DIR}/include/
           DESTINATION ${CMAKE_INSTALL_PREFIX}/include/modules/
           FILES_MATCHING PATTERN "*d")
endif()

# Documentation
option(ENABLE_GREENX_DOCS "Enable documentation" OFF)
if (${ENABLE_GREENX_DOCS})
  find_program(DOXYGEN doxygen REQUIRED)
  message("-- Doxygen documentation support enabled.")
  set(doxy_output "documentation")
  file(MAKE_DIRECTORY ${doxy_output})
  add_custom_target(docs COMMAND ${DOXYGEN} ${PROJECT_SOURCE_DIR}/Doxyfile)
endif ()

# Set installation location
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(default_install_prefix "${PROJECT_SOURCE_DIR}/install")
  set(CMAKE_INSTALL_PREFIX ${default_install_prefix}
    CACHE STRING "Choose the installation directory. Default location is ${default_install_prefix}"
    FORCE)
endif()

install(
    EXPORT      greenXTargets
    NAMESPACE   greenX::
    FILE        greenXTargets.cmake
    DESTINATION lib/cmake/greenX
)

install(
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/greenXConfig.cmake"
    DESTINATION lib/cmake/greenX
)
