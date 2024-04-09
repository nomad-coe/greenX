cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

if(NOT GPU AND NOT NVIDIA)
    project(libpaw_greenx Fortran CXX C)
else()
    project(libpaw_greenx Fortran CXX C)
endif()

# Set the compilers to GNU
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_Fortran_COMPILER "gfortran")

# Set the properties
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -cpp -march=native -fallow-argument-mismatch -ffree-line-length-none -fPIC")

# Where to save the modules
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${MOD_DIR})
if(CMAKE_Fortran_COMPILER MATCHES "ifx" OR CMAKE_Fortran_COMPILER MATCHES "ifort")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -module ${MOD_DIR}")
else()
    set(CMAKE_Fortran_MODULE_DIRECTORY ${MOD_DIR})
endif()

# Some definitions for the compilation to success
add_compile_definitions(HAVE_FC_ISO_C_BINDING)
add_compile_definitions(HAVE_STDDEF_H)
add_compile_definitions(HAVE_STDLIB_H)
add_compile_definitions(HAVE_STDIO_H)
add_compile_definitions(HAVE_MALLOC_H)
add_compile_definitions(HAVE_MATH_H)

# Abinit common files
set(ABINIT_COMMON_SRC "")
add_subdirectory(common/src/02_clib                 )
add_subdirectory(common/src/10_defs                 )
add_subdirectory(common/src/11_memory_mpi           )
add_subdirectory(common/src/12_hide_mpi             )
add_subdirectory(common/src/14_hidewrite            )
add_subdirectory(common/src/16_hideleave            )
add_subdirectory(common/src/17_gpu_toolbox          )
add_subdirectory(common/src/17_libtetra_ext         )
add_subdirectory(common/src/17_minimax              )
add_subdirectory(common/src/17_yaml_out             )
add_subdirectory(common/src/18_timing               )
add_subdirectory(common/src/21_hashfuncs            )
add_subdirectory(common/src/27_toolbox_oop          )
add_subdirectory(common/src/28_numeric_noabirule    )
add_subdirectory(common/src/29_kpoints              )
add_subdirectory(common/src/32_util                 )
add_subdirectory(common/src/33_xc_lowlevel          )

add_library(abinit_common_objects OBJECT ${ABINIT_COMMON_SRC})
target_include_directories(abinit_common_objects PUBLIC common/src/incs)
add_library(abinit_common STATIC $<TARGET_OBJECTS:abinit_common_objects>)

# Now for the libpaw we compile here
add_subdirectory(libpaw/src/)


# Install
set(INSTALL_DIR "${PROJECT_SOURCE_DIR}/libpaw-greenx")
message(STATUS "INSTALL DIRECTORY : ${INSTALL_DIR}")
file(MAKE_DIRECTORY ${INSTALL_DIR})
file(MAKE_DIRECTORY ${INSTALL_DIR}/lib)
# Here we install the modules
install(DIRECTORY ${MOD_DIR}
        DESTINATION ${INSTALL_DIR}
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE
                         GROUP_READ
                         WORLD_READ
        )

# Here we install the libraries
install(TARGETS abinit_common DESTINATION ${INSTALL_DIR}/lib)
install(TARGETS paw DESTINATION ${INSTALL_DIR}/lib)
