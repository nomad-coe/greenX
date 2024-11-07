# ***************************************************************************************************
#  Copyright (C) 2020-2024 GreenX library
#  This file is distributed under the terms of the APACHE2 License.
#
# ***************************************************************************************************

# -----------------------------------------------
# CMakeLists for NAO Integrals Library
# -----------------------------------------------
add_library(LibGXNAO "")

set_target_properties(LibGXNAO
  PROPERTIES
  VERSION 0.0.1
  SOVERSION 0.0.1
  LIBRARY_OUTPUT_NAME gx_nao
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY})

target_include_directories(LibGXNAO PUBLIC src/)

target_sources(LibGXNAO PRIVATE 
    src/math.f90
    src/wigner.f90 
    src/gaunt.f90 
    src/spline.f90
    src/gauss_quadrature.f90
    src/legendre_polynomial.f90
    src/log_grid.f90
    src/codensity_radial_function.f90
    src/spherical_harmonics.f90
    src/codensity_expansion.f90
)
target_link_libraries(LibGXNAO GXCommon ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# -----------------------------------------------
# Library Installation
# -----------------------------------------------
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(TARGETS LibGXNAO ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}  DESTINATION include)

# -----------------------------------------------
# Application Testing Set-Up
# -----------------------------------------------
# Include cmake custom function
include(../cmake/testFunctions.cmake)

# Set name of test sub-directory in the build directory
set(TEST_TARGET_DIR "nao_integrals")

target_include_directories(LibGXNAO PUBLIC test/)

# -----------------------------------------------
# Application Tests
# -----------------------------------------------

# Add a test target
add_executable(test_nao)

# Set binary name
set_target_properties(test_nao
        PROPERTIES
        RUNTIME_OUTPUT_NAME test_nao.exe)

# Define source that comprise the binary
target_sources(test_nao
        PRIVATE
        test/test_nao.f90
        )

# Libraries that the binary links to
target_link_libraries(test_nao
        PUBLIC
	    LibGXNAO
        )

# Build location of the binary
# TODO(Alex) Consider move this to test/, to sit with the python drivers
set_target_properties(test_nao
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})