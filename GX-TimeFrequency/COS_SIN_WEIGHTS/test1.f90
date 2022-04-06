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
      INTEGER, PARAMETER                                 :: unit_nr = 1 
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, PARAMETER                                 :: num_integ_points = 6
      LOGICAL, INTENT(IN)                                :: do_im_time, do_ri_sos_laplace_mp2, &
                                                            do_print
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
!      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_gw_im_time, do_kpoints_cubic_RPA
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: a_scaling
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: e_fermi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: tj, wj
!      TYPE(mp2_type), OPTIONAL, POINTER                  :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_minimax_grid'
      INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200

      INTEGER                                            :: handle, ierr, ispin, jquad, nspins, i_exp!, ii
      LOGICAL                                            :: my_do_kpoints, my_open_shell
      REAL(KIND=dp)                                      :: E_Range, Emax, Emin, max_error_min, &
                                                            scaling, Range_from_i_exp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw


    test1 = get_minimax_grid(homo, Eigenval, &
                               do_im_time, do_ri_sos_laplace_mp2, do_print, tau_tj, tau_wj, do_gw_im_time, &
                               do_kpoints_cubic_RPA, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)
                               

end function 

end module test    

