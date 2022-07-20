cmake_minimum_required(VERSION 3.16)
project(ac_abinit LANGUAGES Fortran VERSION 1.0)

add_library(pade_abinit kinds.f90
            pade_approximant.f90
            test.f90)
add_executable(ac_abinit main.f90)
target_link_libraries(ac_abinit pade_abinit)
