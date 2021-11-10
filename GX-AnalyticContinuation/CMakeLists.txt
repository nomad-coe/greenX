cmake_minimum_required(VERSION 3.16)
project(ac LANGUAGES Fortran VERSION 1.0)

add_library(pade m_simple_pade.f90)
add_executable(ac main.f90)
target_link_libraries(ac pade)
