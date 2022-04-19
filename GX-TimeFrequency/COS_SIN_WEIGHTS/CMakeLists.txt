cmake_minimum_required(VERSION 3.16)
project(ct_st LANGUAGES Fortran VERSION 1.0)

add_library(minimax_ct_st kinds.f90
            mathconstants.f90
            minimax_exp.f90
            minimax_exp_gw.f90
            minimax_rpa.f90
            mp2_grids.f90
            minimax_exp_k15.f90
            minimax_exp_k53.f90)
add_executable(ct_st main.f90)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(ct_st minimax_ct_st  ${BLAS_LIBRARIES}  ${LAPACK_LIBRARIES})
