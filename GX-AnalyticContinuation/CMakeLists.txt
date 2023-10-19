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
target_sources(LibGXAC PRIVATE src/pade_approximant.f90 api/gx_ac.f90)
target_link_libraries(LibGXAC GXCommon)
target_link_libraries(LibGXAC GMP::GMPXX)

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
# Unit Testing Set-up
# -----------------------------------------------
if (${ENABLE_UNITTESTS})
  # Create a directory in the build folder to place generated test drivers and binaries
  set(UNIT_TEST_DIR "${PROJECT_BINARY_DIR}/unit_tests/analytic-continuation")
  file(MAKE_DIRECTORY ${UNIT_TEST_DIR})
  message("-- Analytic continuation unit tests written to: ${UNIT_TEST_DIR}")

  # Libraries on which the unit tests depend
  set(LIBS_FOR_UNIT_TESTS LibGXAC)

  # For unit tests for many modules, put this in a loop
  create_unit_test_executable(TARGET_TEST_DIR ${UNIT_TEST_DIR}
    TEST_NAME "test_pade_approximant"
    REQUIRED_LIBS ${LIBS_FOR_UNIT_TESTS})
endif()