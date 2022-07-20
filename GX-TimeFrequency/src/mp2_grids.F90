!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate frequency and time grids (integration points and weights)
!>        for correlation methods as well as weights for the inhomogeneous cosine/sine transform.
! **************************************************************************************************

#define _REGISTER_EXC(msg) call register_exc(msg, __FILE__, __LINE__)


MODULE mp2_grids

   USE kinds,                           ONLY: dp !, register_exc
   USE constants,                       ONLY: pi
   USE minimax_exp_gw,                  ONLY: get_exp_minimax_coeff_gw
   USE minimax_rpa,                     ONLY: get_rpa_minimax_grids

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_grids'

   PUBLIC :: get_minimax_grid
   PUBLIC :: gx_minimax_grid

CONTAINS

!> \brief Compute minimax grid for RPA energy and GW calculation on imaginary time/frequency domain.
!> \param[in] num_integ_points: Number of mesh points.
!> \param[in] emin: Minimum transition energy.
!> \param[in] emax: Maximum transition energy.
!> \param[out] tau_tj: tau mesh.
!> \param[out] tau_wj: weights for tau mesh.
!> \param[out] tj: imaginary frequency mesh.
!> \param[out] wj: weights for imaginary frequency mesh.
!> \param[out] weights_cos_tf_t_to_w: Weights for tau -> w cosine transform.
!> \param[out] weights_cos_tf_w_to_t: Weights for w -> tau cosine transform.
!> \param[out] weights_sin_tf_t_to_w: Weights for tau -> w sine transform.
!> \param[out] max_error_min: Max error for the three kind of transforms (same order as in args)

subroutine gx_minimax_grid(num_integ_points, emin, emax, &
                           tau_tj, tau_wj, tj, wj, weights_cos_tf_t_to_w, &
                           weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, max_error_min, ierr)

   INTEGER, INTENT(IN)                                :: num_integ_points
   real(dp),intent(in) :: emin, emax
   REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: tau_tj(:), tau_wj(:)
   REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: tj(:), wj(:)
   REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: weights_cos_tf_t_to_w(:,:), &
                                                         weights_cos_tf_w_to_t(:,:), &
                                                         weights_sin_tf_t_to_w(:,:)
   REAL(dp), intent(out) :: max_error_min(3)
   integer, intent(out) :: ierr

   CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gx_minimax_grid'
   INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200

   INTEGER                                            :: jquad
   REAL(dp)                                           :: e_range, scaling
   REAL(dp),allocatable                               :: x_tw(:)

   E_Range = Emax / Emin
   ierr = 0

   ALLOCATE(x_tw(2*num_integ_points))

   CALL get_rpa_minimax_grids(num_integ_points, E_Range, x_tw)
   !if (ierr /= 0) return

   ALLOCATE (tj(num_integ_points))
   ALLOCATE (wj(num_integ_points))

   DO jquad = 1, num_integ_points
      tj(jquad) = x_tw(jquad)
      wj(jquad) = x_tw(jquad + num_integ_points)
   END DO

   ! scale the minimax parameters
   tj(:) = tj(:)*Emin
   wj(:) = wj(:)*Emin

   ! set up the minimax time grid
   CALL get_exp_minimax_coeff_gw(num_integ_points, E_Range, x_tw, ierr)
   if (ierr /= 0) return

   ! For RPA we include already a factor of two (see later steps)
   scaling = 2.0_dp

   ALLOCATE (tau_tj(num_integ_points))
   !ALLOCATE (tau_tj(0:num_integ_points))   CPK original allocation.
   ALLOCATE (tau_wj(num_integ_points))

   DO jquad = 1, num_integ_points
      tau_tj(jquad) = x_tw(jquad)/scaling
      tau_wj(jquad) = x_tw(jquad + num_integ_points)/scaling
   END DO

   ! scale grid from [1,R] to [Emin,Emax]
   tau_tj(:) = tau_tj(:)/Emin
   tau_wj(:) = tau_wj(:)/Emin

   ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
   ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
   ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
   weights_cos_tf_t_to_w = 0.0_dp
   weights_cos_tf_w_to_t = 0.0_dp
   weights_sin_tf_t_to_w = 0.0_dp

   ! get the weights for the cosine transform W^c(it) -> W^c(iw)
   CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                     Emin, Emax, max_error_min(1), num_points_per_magnitude, ierr)
   if (ierr /= 0) return

   ! get the weights for the cosine transform W^c(iw) -> W^c(it)
   CALL get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, tj, &
                                     Emin, Emax, max_error_min(2), num_points_per_magnitude, ierr)
   if (ierr /= 0) return

   ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
   CALL get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, tj, &
                                     Emin, Emax, max_error_min(3), num_points_per_magnitude, ierr)
   if (ierr /= 0) return

   ! Compute the "real" weights used for the inhomogeneous cosine/ FT and check whether
   ! the two matrices for the forward/backward transform are the inverse of each other.
#if 0
   do it=1,gwr%ntau
     do iw=1,gwr%ntau
       wt_cos_wgs(iw, it) = wt_cos_wgs(iw, it) * cos(gwr%tau_mesh(it) * gwr%iw_mesh(iw))
       tw_cos_wgs(it, iw) = tw_cos_wgs(it, iw) * cos(gwr%tau_mesh(it) * gwr%iw_mesh(iw))
     end do
   end do

   ALLOCATE (mat, (gwr%ntau, gwr%ntau))
   mat = matmul(gwr%wt_cos_wgs, gwr%tw_cos_wgs)
   do it=1,gwr%ntau
     mat(it, it) = mat(it, it) - one
   end do
   !print *, "mat", mat
   !cosft_duality_error = maxval(abs(mat))
   DEALLOCATE(mat)
#endif

   DEALLOCATE(x_tw)

end subroutine gx_minimax_grid

! **************************************************************************************************

   SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude, ierr)

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
      integer,intent(out) :: ierr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_t_to_w'

      INTEGER                                            :: iii, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U !, cos_mat, product_mat

      ierr = 0

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
      !ALLOCATE (cos_mat(num_integ_points, num_integ_points))
      !cos_mat = 0.0_dp
      !ALLOCATE (product_mat(num_integ_points, num_integ_points))
      !product_mat = 0.0_dp
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

         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
           _REGISTER_EXC("DGESDD returned ierr != 0")
           return
         end if

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

   END SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w

   SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude, ierr)

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
      integer,intent(out)                                :: ierr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_sin_tf_t_to_w'

      INTEGER                                            :: iii, jjj, jquad, lwork, num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U


      ierr = 0

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
            !y_values(iii) = 2.0_dp*x_values(iii)/((x_values(iii))**2+omega**2)
            y_values(iii) = 2.0_dp*omega/((x_values(iii))**2 + omega**2)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               mat_A(iii, jjj) = SIN(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
            END DO
         END DO

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
           _REGISTER_EXC("DGESDD returned ierr != 0")
           return
         end if

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

   END SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w
!*********************************************************************************************************

   SUBROUTINE get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude, ierr)

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
      integer,intent(out)                                :: ierr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_w_to_t'

      INTEGER                                            :: iii, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega, &
                                                            tau, x_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: omega_wj_work, sing_values, vec_UTy, &
                                                            work, work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      ierr = 0

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

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
           _REGISTER_EXC("DGESDD returned ierr != 0")
           return
         end if

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
      INTENT(IN)                                         :: omega_tj, omega_wj_work, x_values, &
                                                            y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_omega_grid_with_cosine'

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

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

   END SUBROUTINE calc_max_error_fit_omega_grid_with_cosine

   SUBROUTINE get_minimax_grid(num_integ_points, &
                               tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w,ounit)

      INTEGER, PARAMETER                                 :: unit_nr = 1
      INTEGER, INTENT(IN), OPTIONAL                      :: ounit
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: a_scaling
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: e_fermi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_minimax_grid'
      INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200

      INTEGER                                            :: ierr, jquad, i_exp
      REAL(KIND=dp)                                      :: E_Range, Emax, Emin, max_error_min, &
                                                            scaling, Range_from_i_exp

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw



        ierr = 0
        ALLOCATE (tj(num_integ_points))
        ALLOCATE (wj(num_integ_points))
        ALLOCATE (tau_tj(0:num_integ_points))
        ALLOCATE (tau_wj(num_integ_points))
        ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
        ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
        ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
        ALLOCATE (x_tw(2*num_integ_points))
        x_tw = 0.0_dp
         ! @Maryam: Start of print statements that have been added

         DO i_exp = 0, 50

           Range_from_i_exp = 1.58489319246_dp**i_exp

           !IF (unit_nr > 0) THEN
             IF (present(ounit)) then
              ! GET AND PRINT FREQUENCY GRIDS

              IF (num_integ_points .LE. 5) THEN
                 write(*,*)"The grid size you choose is not available."
              ELSE
                 CALL get_rpa_minimax_grids(num_integ_points, Range_from_i_exp, x_tw)
              END IF

              WRITE (ounit, FMT="(T3,A,T66,F15.4)") "Range for the minimax approximation:", Range_from_i_exp
              WRITE (ounit, FMT="(T3,A)") "minimax frequency omega_i     weight of frequency omega_i"
              DO jquad = 1, num_integ_points
                 WRITE (ounit, FMT="(T15,F20.10,F20.10)") x_tw(jquad), x_tw(jquad + num_integ_points)
              END DO

              ! GET AND PRINT TIME GRIDS

              IF (num_integ_points .LE. 5) THEN
                 write(*,*)"The grid size you choose is not available."
              ELSE
                 CALL get_exp_minimax_coeff_gw(num_integ_points, Range_from_i_exp, x_tw, ierr)
              END IF

              WRITE (ounit, FMT="(T3,A)") "minimax time tau_j     weight of time tau_j"
              DO jquad = 1, num_integ_points
                 WRITE (ounit, FMT="(T15,F20.10,F20.10)") x_tw(jquad), x_tw(jquad + num_integ_points)
              END DO

              WRITE (ounit, FMT="(T3,A)") " "


           END IF

         END DO

         ! @Maryam: End of print statements that have been added

         DEALLOCATE (x_tw)

       DO i_exp = 0, 50

          Range_from_i_exp = 1.58489319246_dp**i_exp

          Emin = 1.0_dp
          Emax = Range_from_i_exp
          E_Range = Emax/Emin



!         IF (.NOT. do_ri_sos_laplace_mp2) THEN
             ALLOCATE (x_tw(2*num_integ_points))
             x_tw = 0.0_dp
             tj = 0.0_dp
             wj = 0.0_dp
             tau_tj = 0.0_dp
             tau_wj = 0.0_dp
             weights_cos_tf_t_to_w = 0.0_dp
             weights_cos_tf_w_to_t = 0.0_dp
             weights_sin_tf_t_to_w = 0.0_dp

             ierr = 0
             IF (num_integ_points .LE. 5) THEN
                write(*,*)"The grid size you choose is not available."
             ELSE
                CALL get_rpa_minimax_grids(num_integ_points, E_Range, x_tw)
             END IF



             DO jquad = 1, num_integ_points
                tj(jquad) = x_tw(jquad)
                wj(jquad) = x_tw(jquad + num_integ_points)
             END DO

             DEALLOCATE (x_tw)



             ! scale the minimax parameters
             tj(:) = tj(:)*Emin
             wj(:) = wj(:)*Emin


      ! set up the minimax time grid
 !     IF (do_im_time .OR. do_ri_sos_laplace_mp2) THEN
         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp


         IF (num_integ_points .LE. 5) THEN
            write(*,*)"The grid size you choose is not available."
         ELSE
            CALL get_exp_minimax_coeff_gw(num_integ_points, E_Range, x_tw, ierr)
         END IF

         ! For RPA we include already a factor of two (see later steps)
         scaling = 2.0_dp
!         IF (do_ri_sos_laplace_mp2) scaling = 1.0_dp



         DO jquad = 1, num_integ_points
            tau_tj(jquad) = x_tw(jquad)/scaling
            tau_wj(jquad) = x_tw(jquad + num_integ_points)/scaling
         END DO

         DEALLOCATE (x_tw)

         ! scale grid from [1,R] to [Emin,Emax]
         !tau_tj(:) = tau_tj(:)/Emin*Emin
         !tau_wj(:) = tau_wj(:)/Emin*Emin

         CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                           Emin, Emax, max_error_min, num_points_per_magnitude, ierr)

         CALL get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, tj, &
                                           Emin, Emax, max_error_min, num_points_per_magnitude, ierr)

             !write(1,*) "I am writing tj", tj, "***************************"


         CALL get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, tj, &
                                           Emin, Emax, max_error_min, num_points_per_magnitude, ierr)

             IF (present(ounit)) then
                write(ounit, FMT="(T3,A,T66,F15.4)") "Range for the minimax approximation:", Range_from_i_exp
                write(ounit,*)Range_from_i_exp, "weights_cos_tf_t_to_w", weights_cos_tf_t_to_w
                write(ounit,*)Range_from_i_exp, "weights_cos_tf_w_to_t", weights_cos_tf_w_to_t
                write(ounit,*)Range_from_i_exp, "weights_sin_tf_t_to_w", weights_sin_tf_t_to_w
             ENDIF
                 !IF (unit_nr > 0) THEN
                     !WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
                    ! "MINIMAX_INFO| Maximum deviation of the imag. time fit:", max_error_min
                 !END IF

             END DO

         !DEALLOCATE(tj, wj, tau_tj, tau_wj, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)


   END SUBROUTINE get_minimax_grid

END MODULE mp2_grids
