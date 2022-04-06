module test

  use mp2_grids
  
   USE kinds,                           ONLY: dp
!   USE kpoint_types,                    ONLY: get_kpoint_info,&
!                                              kpoint_env_type,&
!                                              kpoint_type
!   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
!   USE message_passing,                 ONLY: mp_bcast,&
!                                              mp_max,&
!                                              mp_sum
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
   USE minimax_exp_gw,                  ONLY: get_exp_minimax_coeff_gw
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff,&
                                              get_rpa_minimax_coeff_larger_grid  

  implicit none

  public :: test1

contains

function test1()

    complex:: test1
      INTEGER, PARAMETER                                 :: num_integ_points = 6

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj
!      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tj, wj
!      TYPE(mp2_type), OPTIONAL, POINTER                  :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_minimax_grid'
      INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200
      REAL(KIND=dp)                                      :: a_scaling
      INTEGER                                            :: handle, ierr, ispin, jquad, nspins, i_exp!, ii
      REAL(KIND=dp)                                      :: E_Range, Emax, Emin, max_error_min, &
                                                            scaling, Range_from_i_exp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw
      REAL(KIND=dp)                                      :: e_fermi

    call get_minimax_grid(num_integ_points, &
                               tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)
!    test1 = get_minimax_grid(num_integ_points, &
!                               tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
!                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)
                               

end function 

end module test    

