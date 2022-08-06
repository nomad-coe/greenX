# CMakeLists for TimeFrequency Library (MiniMax)
# TODO(Alex)This file should have the option to build/install the minimax library

# GreenXMiniMax Library
add_library(LibGXMiniMax "")

# Set properties of the library
set_target_properties(LibGXMiniMax
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME gx_minimax
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Include src as a directory visible to CMake
target_include_directories(LibGXMiniMax PUBLIC src/ api/)

# Define source that comprise LibGXMiniMax
target_sources(LibGXMiniMax PRIVATE
        src/gx_common.h
        src/minimax_gw.F90
        src/minimax_rpa.F90
        src/mp2_grids.F90
        api/api_utilities.f90
        api/gx_minimax.f90
        )

# Link external libraries and lower-level GX lib dependency
target_link_libraries(LibGXMiniMax GXCommon LAPACK_LIBRARIES BLAS_LIBRARIES)

# -----------------------------------------------
# Library Installation
# -----------------------------------------------

# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}
install(TARGETS LibGXMiniMax
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION include)

# -----------------------------------------------
# Application Tests
# Note, for the addition of mutliple tests, these
# commands can be abstracted into a function.
# -----------------------------------------------
target_include_directories(LibGXMiniMax PUBLIC test/)

add_executable(test_minimax)

# Set binary name
set_target_properties(test_minimax
        PROPERTIES
        RUNTIME_OUTPUT_NAME test_minimax.exe)

# Define source that comprise the binary
target_sources(test_minimax
        PRIVATE
        test/test_gx_minimax_grid.f90
        )

# Libraries that the binary links to
target_link_libraries(test_minimax
        PUBLIC
        LibGXMiniMax
        LAPACK_LIBRARIES
        BLAS_LIBRARIES
        )

# Build location of the binary
# TODO(Alex) Consider move this to test/, to sit with the python drivers
set_target_properties(test_minimax
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_Fortran_BIN_DIRECTORY})

# Copy .py test to the `build/test` directory, such that one can run pytest there
# where CMAKE_CURRENT_SOURCE_DIR => CMakeLists.txt on this level
add_custom_command(
        TARGET LibGXMiniMax POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                # Test source relative to the time-frequency folder
                ${CMAKE_CURRENT_SOURCE_DIR}/test/test_gx_minimax_grid.py
                # Location to copy the test to
                ${PROJECT_BINARY_DIR}/test/time-frequency/test_gx_minimax_grid.py)
# Note.
# The command below would build the test in <BUILD_DIR>/<LIB_NAME>/test:
#                #  Location to copy the test to
#                ${CMAKE_CURRENT_BINARY_DIR}/test/test_gx_minimax_grid.py)
# whereas I am currently opting to put tests in <BUILD_DIR>/tests/<LIB_NAME>

# Add test to ctest
add_test(
        NAME test_minimax
        COMMAND pytest -s
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/time-frequency
)
