      program  main

         use kinds,          ONLY : dp
         use mp2_grids,      only : get_minimax_grid

      INTEGER, PARAMETER                                 :: num_integ_points = 22
      INTEGER                                            :: ct_unit, ct_main_unit
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w

!      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_minimax_grid'
!      INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200
      REAL(KIND=dp)                                      :: a_scaling
!      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw
      REAL(KIND=dp)                                      :: e_fermi

!             ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
!             weights_cos_tf_t_to_w = 0.0_dp
!             ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
!             weights_cos_tf_w_to_t = 0.0_dp
!             ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
!             weights_sin_tf_t_to_w = 0.0_dp
!             ALLOCATE (x_tw(2*num_integ_points))
!             x_tw = 0.0_dp


    open(newunit=ct_unit,file='CT')
    open(newunit=ct_main_unit,file='CT_main')

!    DO i_exp = 0, 50

!           Range_from_i_exp = 1.58489319246_dp**i_exp

    call get_minimax_grid(num_integ_points, &
                               tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, ounit=ct_unit)


!    write(ct_main_unit,*) weights_cos_tf_w_to_t

!    DEALLOCATE(tj, wj, tau_tj, tau_wj, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)

!    END DO

    close(ct_unit)
    close(ct_main_unit)
      end program
