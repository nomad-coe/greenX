# CMakeLists.txt for GreenX Common Modules, required by all other libs in
# the package


add_library(GXCommon "")

set_target_properties(GXCommon
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME GXCommon
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

#target_include_directories(GXCommon PUBLIC src/)
target_include_directories(GXCommon
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>  # used when building in the source tree
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_Fortran_MODULES}/common/modules>                        # used after 'make install'
)

target_sources(GXCommon PRIVATE
        src/kinds.f90
        src/constants.f90
        src/unit_conversion.f90
        src/error_handling.f90
        src/lapack_interfaces.f90
        )

# Install library
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in the top-level
# CMakeLists.txt
install(TARGETS GXCommon
    EXPORT greenXTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

# Install modules
# Destination relative to ${CMAKE_INSTALL_PREFIX}, defined in the top-level
## CMakeLists.txt
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
        DESTINATION ${CMAKE_INSTALL_Fortran_MODULES}/common)


# ensure that zofu is build before GXCommon 
if (${ZOFU_BUILD_ON_FLY})
  add_dependencies(GXCommon zofu)
endif()
