      program  main

         use mp2_grids
         use test

      INTEGER, PARAMETER                                 :: num_integ_points = 14

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
!      REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)          :: Eval
!      REAL (KIND=dp)                                     :: homo
      
!      ALLOCATE(Eval(10))
!      Eval = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 /) 
      

!             ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
!             weights_cos_tf_t_to_w = 0.0_dp
!             ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
!             weights_cos_tf_w_to_t = 0.0_dp
!             ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
!             weights_sin_tf_t_to_w = 0.0_dp
!             ALLOCATE (x_tw(2*num_integ_points))
!             x_tw = 0.0_dp             


    open(1,file='CT')

    call get_minimax_grid(num_integ_points, &
                               tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)
    write(*,*) weights_cos_tf_w_to_t

    close(1)
!    DEALLOCATE(Eval)
      end program


    

