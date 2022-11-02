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
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

target_include_directories(LibGXAC PUBLIC src/)

target_sources(LibGXAC PRIVATE
        src/pade_approximant.f90
        )

target_link_libraries(LibGXAC GXCommon)

# -----------------------------------------------
# Library Installation
# -----------------------------------------------
# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(TARGETS LibGXAC
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in top-level CMake
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION include)

# -----------------------------------------------
# Unit Testing Set-up
# -----------------------------------------------

