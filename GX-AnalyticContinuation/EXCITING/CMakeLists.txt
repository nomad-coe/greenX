cmake_minimum_required(VERSION 3.16)
project(ac LANGUAGES Fortran VERSION 1.0)

add_library(pade_exciting constants.f90
            pade.F90
            ctdfrac.F90
            mod_kpoint.F90
            mod_pade.f90
            trace.F90
            precision.f90
            modmpi.F90)
add_executable(ac_exciting main.f90)
target_link_libraries(ac_exciting pade_exciting)
