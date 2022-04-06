MODULE mp2_grids
!   USE cp_para_types,                   ONLY: cp_para_env_type
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
!   USE mp2_types,                       ONLY: mp2_type
!   USE qs_environment_types,            ONLY: get_qs_env,&
!                                              qs_environment_type
!   USE qs_mo_types,                     ONLY: get_mo_set,&
!                                              mo_set_type
   
   IMPLICIT NONE


   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_grids'

   PUBLIC :: get_minimax_grid

CONTAINS


   SUBROUTINE get_minimax_grid(num_integ_points, &
                               tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)



!      TYPE(cp_para_env_type), POINTER                    :: para_env
!      INTEGER, INTENT(IN)                                :: unit_nr = 1
      INTEGER, PARAMETER                                 :: unit_nr = 1 
!      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
!      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points
!      LOGICAL, INTENT(IN)                                :: do_im_time, do_ri_sos_laplace_mp2, &
!                                                            do_print
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
!      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
!      LOGICAL, INTENT(IN), OPTIONAL                      :: do_gw_im_time, do_kpoints_cubic_RPA
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
!      LOGICAL                                            :: my_do_kpoints, my_open_shell
      REAL(KIND=dp)                                      :: E_Range, Emax, Emin, max_error_min, &
                                                            scaling, Range_from_i_exp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw


      ! Test for spin unrestricted
!      nspins = SIZE(homo)
!      my_open_shell = (nspins == 2)

      ! Test whether all necessary variables are available
!      my_do_kpoints = .FALSE.
!      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         !IF (.NOT. (PRESENT(do_gw_im_time) .AND. PRESENT(do_kpoints_cubic_RPA) .AND. PRESENT(a_scaling) &
          !          .AND. PRESENT(e_fermi) .AND. PRESENT(tj) .AND. PRESENT(wj) .AND. PRESENT(weights_cos_tf_t_to_w) .AND. &
          !          PRESENT(weights_cos_tf_w_to_t) .AND. PRESENT(weights_sin_tf_t_to_w) .AND. PRESENT(mp2_env))) THEN
          !  CPABORT("Need more parameters without SOS-MP2")
         !END IF
!         my_do_kpoints = do_kpoints_cubic_RPA
!      END IF


      !Emin = 2.0
      !Emax = 2.0*E_Range
      IF (num_integ_points > 20 .AND. E_Range < 100.0_dp) THEN
         IF (unit_nr > 0) &
           write(1,*)    "You requested a large minimax grid (> 20 points) for a small minimax range R (R < 100). "// &
                         "That may lead to numerical "// &
                         "instabilities when computing minimax grid weights. You can prevent small ranges by choosing "// &
                         "a larger basis set with higher angular momenta or alternatively using all-electron calculations."
      END IF

!      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ierr = 0
         
         ! @Maryam: Start of print statements that have been added

         DO i_exp = 0, 50

           Range_from_i_exp = 1.58489319246_dp**i_exp

           IF (unit_nr > 0) THEN

              ! GET AND PRINT FREQUENCY GRIDS

              IF (num_integ_points .LE. 20) THEN
                 CALL get_rpa_minimax_coeff(num_integ_points, Range_from_i_exp, x_tw, ierr)
              ELSE
                 CALL get_rpa_minimax_coeff_larger_grid(num_integ_points, Range_from_i_exp, x_tw)
              END IF

              WRITE (1, FMT="(T3,A,T66,F15.4)") "Range for the minimax approximation:", Range_from_i_exp
              WRITE (1, FMT="(T3,A)") "minimax frequency omega_i     weight of frequency omega_i"
              DO jquad = 1, num_integ_points
                 WRITE (1, FMT="(T15,F20.10,F20.10)") x_tw(jquad), x_tw(jquad + num_integ_points)
              END DO

              ! GET AND PRINT TIME GRIDS

              IF (num_integ_points .LE. 20) THEN
                 CALL get_exp_minimax_coeff(num_integ_points, Range_from_i_exp, x_tw)
              ELSE
                 CALL get_exp_minimax_coeff_gw(num_integ_points, Range_from_i_exp, x_tw)
              END IF

              WRITE (1, FMT="(T3,A)") "minimax time tau_j     weight of time tau_j"
              DO jquad = 1, num_integ_points
                 WRITE (1, FMT="(T15,F20.10,F20.10)") x_tw(jquad), x_tw(jquad + num_integ_points)
              END DO

              WRITE (1, FMT="(T3,A)") " "


           END IF

         END DO

         ! @Maryam: End of print statements that have been added
         
         
         IF (num_integ_points .LE. 20) THEN
            CALL get_rpa_minimax_coeff(num_integ_points, E_Range, x_tw, ierr)
         ELSE
            CALL get_rpa_minimax_coeff_larger_grid(num_integ_points, E_Range, x_tw)
         END IF

         ALLOCATE (tj(num_integ_points))
         tj = 0.0_dp

         ALLOCATE (wj(num_integ_points))
         wj = 0.0_dp

         DO jquad = 1, num_integ_points
            tj(jquad) = x_tw(jquad)
            wj(jquad) = x_tw(jquad + num_integ_points)
         END DO

         DEALLOCATE (x_tw)

!         IF (unit_nr > 0 .AND. do_print) THEN
            WRITE (1, FMT="(T3,A,T75,i6)") &
               "MINIMAX_INFO| Number of integration points:", num_integ_points
            WRITE (1, FMT="(T3,A,T66,F15.4)") &
               "MINIMAX_INFO| Range for the minimax approximation:", E_Range
            WRITE (1, FMT="(T3,A,T54,A,T72,A)") "MINIMAX_INFO| Minimax parameters:", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (1, FMT="(T41,F20.10,F20.10)") wj(jquad), tj(jquad)
            END DO
!            CALL m_flush(unit_nr)
!         END IF

         ! scale the minimax parameters
         tj(:) = tj(:)*Emin
         wj(:) = wj(:)*Emin
!      ELSE
         ! When we perform SOS-MP2, we need an additional factor of 2 for the energies (compare with mp2_laplace.F)
         ! We do not need weights etc. for the cosine transform
         ! We do not scale Emax because it is not needed for SOS-MP2
         Emin = Emin*2.0_dp
!      END IF

      ! set up the minimax time grid
!      IF (do_im_time .OR. do_ri_sos_laplace_mp2) THEN

         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp

         IF (num_integ_points .LE. 20) THEN
            CALL get_exp_minimax_coeff(num_integ_points, E_Range, x_tw)
         ELSE
            CALL get_exp_minimax_coeff_gw(num_integ_points, E_Range, x_tw)
         END IF

         ! For RPA we include already a factor of two (see later steps)
         scaling = 2.0_dp
!         IF (do_ri_sos_laplace_mp2) scaling = 1.0_dp

         ALLOCATE (tau_tj(0:num_integ_points))
         tau_tj = 0.0_dp

         ALLOCATE (tau_wj(num_integ_points))
         tau_wj = 0.0_dp

         DO jquad = 1, num_integ_points
            tau_tj(jquad) = x_tw(jquad)/scaling
            tau_wj(jquad) = x_tw(jquad + num_integ_points)/scaling
         END DO

         DEALLOCATE (x_tw)

!         IF (unit_nr > 0 .AND. do_print) THEN
            WRITE (1, FMT="(T3,A,T66,F15.4)") &
               "MINIMAX_INFO| Range for the minimax approximation:", E_Range
            ! For testing the gap
            WRITE (1, FMT="(T3,A,T66,F15.4)") &
               "MINIMAX_INFO| Gap:", Emin
            WRITE (1, FMT="(T3,A,T54,A,T72,A)") &
               "MINIMAX_INFO| Minimax parameters of the time grid:", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (1, FMT="(T41,F20.10,F20.10)") tau_wj(jquad), tau_tj(jquad)
            END DO
!            CALL m_flush(unit_nr)
!         END IF

         ! scale grid from [1,R] to [Emin,Emax]
         !tau_tj(:) = tau_tj(:)/Emin
         !tau_wj(:) = tau_wj(:)/Emin

!         IF (.NOT. do_ri_sos_laplace_mp2) THEN
             ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
             weights_cos_tf_t_to_w = 0.0_dp
             ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
             weights_cos_tf_w_to_t = 0.0_dp
             ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
             weights_sin_tf_t_to_w = 0.0_dp
             DO i_exp = 0, 50
                 Range_from_i_exp = 1.58489319246_dp**i_exp
                 Emin = 1.0
                 Emax = Range_from_i_exp

                 CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                              Emin, Emax, max_error_min, num_points_per_magnitude)
!                 IF (do_gw_im_time) THEN




                     CALL get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, tj, &
                                                 Emin, Emax, max_error_min, num_points_per_magnitude)


                     write(1,*) "*************************", weights_cos_tf_w_to_t, "***************************"
                     
                     
                     CALL get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, tj, &
                                                 Emin, Emax, max_error_min, num_points_per_magnitude)
                                                 
                                                 
                                                 
                     !WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") "Range for the minimax approximation:",&
                                            !  Range_from_i_exp
                     !WRITE (UNIT=unit_nr, FMT="(T3,A)") " weights"


                     !DO jquad = 1, num_integ_points

                         !WRITE (UNIT=unit_nr, FMT="(T15,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)") Emin, Emax,&
                             !Emax/Emin, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w


                     !END DO

!                 END IF
                 IF (unit_nr > 0) THEN
                     WRITE (1, FMT="(T3,A,T66,ES15.2)") &
                     "MINIMAX_INFO| Maximum deviation of the imag. time fit:", max_error_min
                 END IF


             END DO


!         END IF

!      END IF

   END SUBROUTINE get_minimax_grid
   
   
! **************************************************************************************************


   
   SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_t_to_w'

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U, cos_mat

      

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (cos_mat(num_integ_points, num_integ_points))
      cos_mat = 0.0_dp      
      ALLOCATE (tau_wj_work(num_integ_points))
      tau_wj_work = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      max_error = 0.0_dp

      ! loop over all omega frequency points
      DO jquad = 1, num_integ_points

         ! set the x-values logarithmically in the interval [Emin,Emax]
         multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
         DO iii = 1, num_x_nodes
            x_values(iii) = E_min*multiplicator**(iii - 1)
         END DO

         omega = omega_tj(jquad)

         ! y=2x/(x^2+omega_k^2)
         DO iii = 1, num_x_nodes
            y_values(iii) = 2.0_dp*x_values(iii)/((x_values(iii))**2 + omega**2)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               mat_A(iii, jjj) = COS(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
            END DO
         END DO
         
         DO jjj = 1, num_integ_points
            cos_mat(jquad, jjj) = COS(omega*tau_tj(jjj))
         END DO
         
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

!         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, tau_wj_work, num_integ_points)

         weights_cos_tf_t_to_w(jquad, :) = tau_wj_work(:)

         CALL calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                      y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

      write(1,*) "Emax ", E_max,  "Emin", E_min, "Range for weights", E_max/E_min
      write(1,*) "cos_mat", cos_mat
      write(1,*) 'weights_cos_tf_t_to_w =', weights_cos_tf_t_to_w
      DEALLOCATE (cos_mat)

   END SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w

   SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_sin_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(OUT)                         :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_sin_tf_t_to_w'

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U, sin_mat

      

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (sin_mat(num_integ_points, num_integ_points))
      sin_mat = 0.0_dp      
      ALLOCATE (tau_wj_work(num_integ_points))
      tau_wj_work = 0.0_dp
      ALLOCATE (work_array(2*num_integ_points))
      work_array = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      max_error = 0.0_dp

      ! loop over all omega frequency points
      DO jquad = 1, num_integ_points

         chi2_min_jquad = 100.0_dp

         ! set the x-values logarithmically in the interval [Emin,Emax]
         multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
         DO iii = 1, num_x_nodes
            x_values(iii) = E_min*multiplicator**(iii - 1)
         END DO

         omega = omega_tj(jquad)

         ! y=2x/(x^2+omega_k^2)
         DO iii = 1, num_x_nodes
!            y_values(iii) = 2.0_dp*x_values(iii)/((x_values(iii))**2+omega**2)
            y_values(iii) = 2.0_dp*omega/((x_values(iii))**2 + omega**2)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               mat_A(iii, jjj) = SIN(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
            END DO
         END DO
         
         DO jjj = 1, num_integ_points
            sin_mat(jquad, jjj) = SIN(omega*tau_tj(jjj))
         END DO   
         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

!         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, tau_wj_work, num_integ_points)

         weights_sin_tf_t_to_w(jquad, :) = tau_wj_work(:)

         CALL calc_max_error_fit_tau_grid_with_sine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                    y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, work_array, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)
      write(1,*) "Emax ", E_max,  "Emin", E_min, "Range for weights", E_max/E_min
      write(1,*) "sin_mat" , sin_mat
      write(1,*) 'weights_sin_tf_t_to_w =', weights_sin_tf_t_to_w
      DEALLOCATE (sin_mat)

   END SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w
!*********************************************************************************************************
   SUBROUTINE get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_w_to_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_w_to_t'

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega, &
                                                            tau, x_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: omega_wj_work, sing_values, vec_UTy, &
                                                            work, work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U, cos_mat_w_to_t

      

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (cos_mat_w_to_t(num_integ_points, num_integ_points))
      cos_mat_w_to_t = 0.0_dp
      ALLOCATE (omega_wj_work(num_integ_points))
      omega_wj_work = 0.0_dp
      ALLOCATE (work_array(2*num_integ_points))
      work_array = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      ! set the x-values logarithmically in the interval [Emin,Emax]
      multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
      DO iii = 1, num_x_nodes
         x_values(iii) = E_min*multiplicator**(iii - 1)
      END DO

      max_error = 0.0_dp

      ! loop over all tau time points
      DO jquad = 1, num_integ_points

         chi2_min_jquad = 100.0_dp

         tau = tau_tj(jquad)

         ! y=exp(-x*|tau_k|)
         DO iii = 1, num_x_nodes
            y_values(iii) = EXP(-x_values(iii)*tau)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               omega = omega_tj(jjj)
               x_value = x_values(iii)
               mat_A(iii, jjj) = COS(tau*omega)*2.0_dp*x_value/(x_value**2 + omega**2)
            END DO
         END DO
         DO jjj = 1, num_integ_points
               omega = omega_tj(jjj)
               cos_mat_w_to_t(jquad, jjj) = COS(tau*omega)
         END DO

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

!         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, omega_wj_work, num_integ_points)

         weights_cos_tf_w_to_t(jquad, :) = omega_wj_work(:)

         CALL calc_max_error_fit_omega_grid_with_cosine(max_error, tau, omega_tj, omega_wj_work, x_values, &
                                                        y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, omega_wj_work, work_array, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)
      write(1,*) "Emax ", E_max,  "Emin", E_min, "Range for weights", E_max/E_min
      write(1,*) "cos_mat_w_to_t" , cos_mat_w_to_t
      write(1,*) 'weights_cos_tf_w_to_t =', weights_cos_tf_w_to_t
      DEALLOCATE (cos_mat_w_to_t)

   END SUBROUTINE get_l_sq_wghts_cos_tf_w_to_t   
   
   
! **************************************************************************************************
   PURE SUBROUTINE calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                           y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_cosine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk) - func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk) - func_val)
            func_val_temp = func_val
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

   END SUBROUTINE calc_max_error_fit_tau_grid_with_cosine   
   
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_tau_grid_cosine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp), INTENT(OUT)                         :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: omega

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val + tau_wj_work(iii)*COS(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE eval_fit_func_tau_grid_cosine

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for sine transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_tau_grid_sine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp), INTENT(INOUT)                       :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(in)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: omega

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val + tau_wj_work(iii)*SIN(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE eval_fit_func_tau_grid_sine

! **************************************************************************************************

   PURE SUBROUTINE calc_max_error_fit_tau_grid_with_sine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                         y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_sine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk) - func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk) - func_val)
            func_val_temp = func_val
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

   END SUBROUTINE calc_max_error_fit_tau_grid_with_sine

! **************************************************************************************************

! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_omega_grid_cosine(func_val, x_value, num_integ_points, omega_tj, omega_wj_work, tau)
      REAL(KIND=dp), INTENT(OUT)                         :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj, omega_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: tau

      INTEGER                                            :: iii
      REAL(KIND=dp)                                      :: omega

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         omega = omega_tj(iii)
         func_val = func_val + omega_wj_work(iii)*COS(tau*omega)*2.0_dp*x_value/(x_value**2 + omega**2)

      END DO

   END SUBROUTINE eval_fit_func_omega_grid_cosine



   SUBROUTINE calc_max_error_fit_omega_grid_with_cosine(max_error, tau, omega_tj, omega_wj_work, x_values, &
                                                        y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      REAL(KIND=dp), INTENT(IN)                          :: tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj, omega_wj_work, x_values, &
                                                            y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_omega_grid_with_cosine'

      INTEGER                                            :: handle, kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

!      CALL timeset(routineN, handle)

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_omega_grid_cosine(func_val, x_values(kkk), num_integ_points, omega_tj, omega_wj_work, tau)

         IF (ABS(y_values(kkk) - func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk) - func_val)
            func_val_temp = func_val
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

!      CALL timestop(handle)

   END SUBROUTINE calc_max_error_fit_omega_grid_with_cosine
END MODULE mp2_grids
