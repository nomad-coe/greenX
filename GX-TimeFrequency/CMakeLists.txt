# CMakeLists for TimeFrequency Libraries 
# TODO(Alex) Have this as a separate CMake 
# This file should just have the option to build/install specific libs 

# Minimax Grids Library
add_library(GXMiniMax "")

set_target_properties(GXMiniMax
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.01
        LIBRARY_OUTPUT_NAME GXMiniMax
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_link_libraries(GXMiniMax GXCommon)

target_include_directories(GXMiniMax PUBLIC src/ )

# Set full path to minimax grids using preprocessor 
set_property(SOURCE src/minimax/minimax.F90 APPEND PROPERTY
    COMPILE_DEFINITIONS "MINIMAX_PATH='${CMAKE_SOURCE_DIR}/GX-TimeFrequency/data/minimax_grids'")

target_sources(GXMiniMax PRIVATE
        src/minimax/minimax.F90
        src/minimax/parse_minimax.f90
        src/minimax/read_cmd_line.f90  # TODO The test depends on this, not the lib
        )


# TODO (Alex) This wants to be defined once
# Set installation location
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # TODO(Alex) clean this up
    set(default_install_prefix "${PROJECT_SOURCE_DIR}/install")        
    set(CMAKE_INSTALL_PREFIX ${default_install_prefix}
            CACHE STRING "Choose the installation directory. Default location is ${default_install_prefix}"
            FORCE)
endif()

# Destination relative to ${CMAKE_INSTALL_PREFIX}
install(TARGETS GXMiniMax
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Destination relative to ${CMAKE_INSTALL_PREFIX}
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION include)

# Application tests
set(REG_TESTS test_minimax)

foreach(INPUT ${REG_TESTS})

    # Fortran test TODO(Alex) Move where this is built
    add_executable(${INPUT}.x src/minimax/${INPUT}.f90)
    link_directories(${LIBRARY_OUTPUT_PATH})
    target_link_libraries(${INPUT}.x GXMiniMax)

    # Python test framework
    add_test(
            NAME REG_TEST_${INPUT}
            COMMAND pytest -s GX-TimeFrequency/tests/test_minimax/${INPUT}.py --root ${CMAKE_SOURCE_DIR} --binary ${CMAKE_BINARY_DIR}/GX-TimeFrequency/${INPUT}.x
            # Run in GreenX's root
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

    # Can also set useful test properties
endforeach()

