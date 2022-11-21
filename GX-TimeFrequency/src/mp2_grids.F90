!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate frequency and time grids (integration points and weights)
!> for correlation methods as well as weights for the inhomogeneous cosine/sine transform.
!>
!> NB: When dealing with runtime exceptions, we set ierr to a non-zero value and return immediately
!  to the caller so we don't need to goto to a cleanup section at the end of the procedure.
!  Assume -std=f2008: i.e. allocatable arrays are automatically deallocated when going out of scope.
! **************************************************************************************************

module mp2_grids
#include "gx_common.h"
   use kinds, only: dp
   use error_handling, only: register_exc
   use constants, only: pi
   use minimax_gw, only: get_exp_minimax_coeff_gw
   use minimax_rpa, only: get_rpa_minimax_grids

   implicit none

   private

   !> Main entry point for client code.
   public :: gx_minimax_grid
   !> TODO(Maryam) Temporary API - to delete
   public :: get_minimax_grid

contains

!> \brief Compute minimax grid for RPA energy and GW calculation on imaginary time/frequency domain.
!> \param[in] num_integ_points: Number of mesh points.
!> \param[in] emin: Minimum transition energy. Arbitrary units as we only need emax/emin
!> \param[in] emax: Maximum transition energy.
!> \param[out] tau_mesh: tau mesh (tau_tj in CP2K)
!> \param[out] tau_wgs: weights for tau mesh (tau_wj in CP2K)
!> \param[out] iw_mesh: imaginary frequency mesh (tj in CP2K)
!> \param[out] iw_wgs: weights for imaginary frequency mesh (wj in CP2K)
!> \param[out] cosft_wt: weights for tau -> w cosine transform. cos(w*t) factor is included.
!> \param[out] cosft_tw: weights for w -> tau cosine transform. cos(w*t) factor is included.
!> \param[out] sinft_wt: weights for tau -> w sine transform. sin(w*t) factor is included.
!> \param[out] max_errors: Max error for the three kind of transforms (same order as previous args)
!> \param[out] cosft_duality_error. Max_{ij} |AB - I| where A and B are the cosft_wt and cosft_tw matrices.
!> \param[out] ierr: Exit status
   subroutine gx_minimax_grid(num_integ_points, emin, emax, &
                              tau_mesh, tau_wgs, iw_mesh, iw_wgs, &
                              cosft_wt, cosft_tw, sinft_wt, &
                              max_errors, cosft_duality_error, ierr)

      integer, intent(in) :: num_integ_points
      real(dp), intent(in) :: emin, emax
      real(dp), allocatable, intent(out) :: tau_mesh(:), tau_wgs(:)
      real(dp), allocatable, intent(out) :: iw_mesh(:), iw_wgs(:)
      real(dp), allocatable, intent(out) :: cosft_wt(:, :), cosft_tw(:, :), sinft_wt(:, :)
      real(dp), intent(out) :: max_errors(3), cosft_duality_error
      integer, intent(out) :: ierr

      integer, parameter                                 :: num_points_per_magnitude = 200
      integer                                            :: jquad, it, iw
      real(dp)                                           :: e_range, scaling
      real(dp), allocatable                              :: x_tw(:), mat(:, :)

      E_Range = Emax/Emin
      ierr = 0

      allocate (x_tw(2*num_integ_points))

      call get_rpa_minimax_grids(num_integ_points, E_Range, x_tw, ierr)
      if (ierr /= 0) return

      allocate (iw_mesh(num_integ_points))
      allocate (iw_wgs(num_integ_points))

      do jquad = 1, num_integ_points
         iw_mesh(jquad) = x_tw(jquad)
         iw_wgs(jquad) = x_tw(jquad + num_integ_points)
      end do

      ! scale the minimax parameters
      iw_mesh(:) = iw_mesh(:)*Emin
      iw_wgs(:) = iw_wgs(:)*Emin

      ! set up the minimax time grid
      call get_exp_minimax_coeff_gw(num_integ_points, E_Range, x_tw, ierr)
      if (ierr /= 0) return

      ! For RPA we include already a factor of two (see later steps)
      scaling = 2.0d0

      allocate (tau_mesh(num_integ_points))
      !ALLOCATE (tau_mesh(0:num_integ_points))   CPK original allocation.
      allocate (tau_wgs(num_integ_points))

      do jquad = 1, num_integ_points
         tau_mesh(jquad) = x_tw(jquad)/scaling
         tau_wgs(jquad) = x_tw(jquad + num_integ_points)/scaling
      end do

      ! scale grid from [1,R] to [Emin,Emax]
      tau_mesh(:) = tau_mesh(:)/Emin
      tau_wgs(:) = tau_wgs(:)/Emin

      allocate (cosft_wt(num_integ_points, num_integ_points))
      allocate (cosft_tw(num_integ_points, num_integ_points))
      allocate (sinft_wt(num_integ_points, num_integ_points))

      ! get the weights for the cosine transform W^c(it) -> W^c(iw)
      call get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_mesh, iw_mesh, cosft_wt, &
                                        Emin, Emax, max_errors(1), num_points_per_magnitude, ierr)
      if (ierr /= 0) return

      ! get the weights for the cosine transform W^c(iw) -> W^c(it)
      call get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_mesh, iw_mesh, cosft_tw, &
                                        Emin, Emax, max_errors(2), num_points_per_magnitude, ierr)
      if (ierr /= 0) return

      ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
      call get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_mesh, iw_mesh, sinft_wt, &
                                        Emin, Emax, max_errors(3), num_points_per_magnitude, ierr)
      if (ierr /= 0) return

      ! Compute the actual weights used for the inhomogeneous cosine/ FT and check whether
      ! the two matrices for the forward/backward transform are the inverse of each other.
      do it = 1, num_integ_points
         do iw = 1, num_integ_points
            cosft_wt(iw, it) = cosft_wt(iw, it)*cos(tau_mesh(it)*iw_mesh(iw))
            cosft_tw(it, iw) = cosft_tw(it, iw)*cos(tau_mesh(it)*iw_mesh(iw))
            sinft_wt(iw, it) = sinft_wt(iw, it)*sin(tau_mesh(it)*iw_mesh(iw))
         end do
      end do

      allocate (mat(num_integ_points, num_integ_points))
      mat(:, :) = matmul(cosft_wt, cosft_tw)
      do it = 1, num_integ_points
         mat(it, it) = mat(it, it) - 1.0d0
      end do
      cosft_duality_error = maxval(abs(mat))

      deallocate (mat)
      deallocate (x_tw)

   end subroutine gx_minimax_grid


   subroutine get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, omega_tj, weights, &
                                           E_min, E_max, max_error, num_points_per_magnitude, ierr)

      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj
       real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: omega_tj
      real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                   :: weights
      real(kind=dp), intent(in)                          :: E_min, E_max
      real(kind=dp), intent(inout)                       :: max_error
      integer, intent(in)                                :: num_points_per_magnitude
      integer, intent(out)                               :: ierr

      integer                                            :: iii, jjj, jquad, lwork, num_x_nodes
      integer, allocatable, dimension(:)                 :: iwork
      real(kind=dp)                                      :: multiplicator, current_point
      real(kind=dp), allocatable, dimension(:)           :: weights_work, sing_values, vec_UTy, work, &
                                                            x_values, y_values
      real(kind=dp), allocatable, dimension(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      ierr = 0

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (int(log10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = max(num_x_nodes, num_integ_points)

      allocate (x_values(num_x_nodes))
      x_values = 0.0d0
      allocate (y_values(num_x_nodes))
      y_values = 0.0d0
      allocate (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0d0
      allocate (weights_work(num_integ_points))
      weights_work = 0.0d0
      allocate (sing_values(num_integ_points))
      sing_values = 0.0d0
      allocate (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0d0
      allocate (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0d0
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      allocate (work(lwork))
      work = 0.0d0
      allocate (iwork(8*num_integ_points))
      iwork = 0
      allocate (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0d0
      allocate (vec_UTy(num_x_nodes))
      vec_UTy = 0.0d0

      ! set the x-values logarithmically in the interval [Emin,Emax]
      multiplicator = (E_max/E_min)**(1.0d0/(real(num_x_nodes, kind=dp) - 1.0d0))
      do iii = 1, num_x_nodes
         x_values(iii) = E_min*multiplicator**(iii - 1)
      end do

      max_error = 0.0d0

      ! loop over all omega frequency points
      do jquad = 1, num_integ_points

         current_point = omega_tj(jquad)

         ! y=2x/(x^2+omega_k^2)
         do iii = 1, num_x_nodes
            y_values(iii) = 2.0d0*x_values(iii)/((x_values(iii))**2 + current_point**2)
         end do

         ! calculate mat_A
         do jjj = 1, num_integ_points
            do iii = 1, num_x_nodes
               mat_A(iii, jjj) = cos(current_point*tau_tj(jjj))*exp(-x_values(iii)*tau_tj(jjj))
            end do
         end do

         ! Singular value decomposition of mat_A
         call dgesdd('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
            _REGISTER_EXC("DGESDD returned ierr != 0")
            return
         end if

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         do jjj = 1, num_integ_points
            do iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            end do
         end do

         ! 2) U^T y
         call dgemm('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0d0, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0d0, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         call dgemm('N', 'N', num_integ_points, 1, num_x_nodes, 1.0d0, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0d0, weights_work, num_integ_points)

         weights(jquad, :) = weights_work(:)

         call calc_max_error_fit_tau_grid_with_cosine(max_error, current_point, tau_tj, weights_work, x_values, &
                                                      y_values, num_integ_points, num_x_nodes)

      end do ! jquad

      deallocate (x_values, y_values, mat_A, weights_work, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

   end subroutine get_l_sq_wghts_cos_tf_t_to_w


   subroutine get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj,omega_tj, weights, &
                                           E_min, E_max, max_error, num_points_per_magnitude, ierr)

      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: omega_tj
       real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                   :: weights
      real(kind=dp), intent(in)                          :: E_min, E_max
      real(kind=dp), intent(out)                         :: max_error
      integer, intent(in)                                :: num_points_per_magnitude
      integer, intent(out)                               :: ierr

      integer                                            :: iii, jjj, jquad, lwork, num_x_nodes
      integer, allocatable, dimension(:)                 :: iwork
      real(kind=dp)                                      :: multiplicator, current_point
      real(kind=dp), allocatable, dimension(:)           :: weights_work, sing_values, vec_UTy, work, &
                                                            x_values, y_values
      real(kind=dp), allocatable, dimension(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      ierr = 0

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (int(log10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = max(num_x_nodes, num_integ_points)

      allocate (x_values(num_x_nodes))
      x_values = 0.0d0
      allocate (y_values(num_x_nodes))
      y_values = 0.0d0
      allocate (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0d0
      allocate (weights_work(num_integ_points))
      weights_work = 0.0d0
      allocate (sing_values(num_integ_points))
      sing_values = 0.0d0
      allocate (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0d0
      allocate (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0d0
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      allocate (work(lwork))
      work = 0.0d0
      allocate (iwork(8*num_integ_points))
      iwork = 0
      allocate (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0d0
      allocate (vec_UTy(num_x_nodes))
      vec_UTy = 0.0d0

      ! set the x-values logarithmically in the interval [Emin,Emax]
      multiplicator = (E_max/E_min)**(1.0d0/(real(num_x_nodes, kind=dp) - 1.0d0))
      do iii = 1, num_x_nodes
         x_values(iii) = E_min*multiplicator**(iii - 1)
      end do 

      max_error = 0.0d0

      ! loop over all omega frequency points
      do jquad = 1, num_integ_points

         current_point = omega_tj(jquad)

         ! y=2*omega_k/(x^2+omega_k^2)
         do iii = 1, num_x_nodes
            y_values(iii) = 2.0d0*current_point/((x_values(iii))**2 + current_point**2)
         end do

         ! calculate mat_A
         do jjj = 1, num_integ_points
            do iii = 1, num_x_nodes
               mat_A(iii, jjj) = sin(current_point*tau_tj(jjj))*exp(-x_values(iii)*tau_tj(jjj))
            end do
         end do

         ! Singular value decomposition of mat_A
         call dgesdd('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
            _REGISTER_EXC("DGESDD returned ierr != 0")
            return
         end if

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         do jjj = 1, num_integ_points
            do iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            end do
         end do

         ! 2) U^T y
         call dgemm('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0d0, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0d0, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         call dgemm('N', 'N', num_integ_points, 1, num_x_nodes, 1.0d0, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0d0, weights_work, num_integ_points)

         weights(jquad, :) = weights_work(:)

         call calc_max_error_fit_tau_grid_with_sine(max_error, current_point, tau_tj, weights_work, x_values, &
                                                    y_values, num_integ_points, num_x_nodes)

      end do ! jquad

      deallocate (x_values, y_values, mat_A, weights_work, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

   end subroutine get_l_sq_wghts_sin_tf_t_to_w


   subroutine get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, omega_tj, weights, &
                                           E_min, E_max, max_error, num_points_per_magnitude, ierr)

      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: omega_tj
       real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                   :: weights
      real(kind=dp), intent(in)                          :: E_min, E_max
      real(kind=dp), intent(inout)                       :: max_error
      integer, intent(in)                                :: num_points_per_magnitude
      integer, intent(out)                               :: ierr

      integer                                            :: iii, jjj, jquad, lwork, num_x_nodes
      integer, allocatable, dimension(:)                 :: iwork
      real(kind=dp)                                      :: multiplicator, current_point, &
                                                            omega,x_value
      real(kind=dp), allocatable, dimension(:)           :: weights_work, sing_values, vec_UTy, &
                                                            work, x_values, y_values
      real(kind=dp), allocatable, dimension(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      ierr = 0

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (int(log10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = max(num_x_nodes, num_integ_points)

      allocate (x_values(num_x_nodes))
      x_values = 0.0d0
      allocate (y_values(num_x_nodes))
      y_values = 0.0d0
      allocate (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0d0
      allocate (weights_work(num_integ_points))
      weights_work = 0.0d0
      allocate (sing_values(num_integ_points))
      sing_values = 0.0d0
      allocate (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0d0
      allocate (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0d0
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      allocate (work(lwork))
      work = 0.0d0
      allocate (iwork(8*num_integ_points))
      iwork = 0
      allocate (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0d0
      allocate (vec_UTy(num_x_nodes))
      vec_UTy = 0.0d0

      ! set the x-values logarithmically in the interval [Emin,Emax]
      multiplicator = (E_max/E_min)**(1.0d0/(real(num_x_nodes, kind=dp) - 1.0d0))
      do iii = 1, num_x_nodes
         x_values(iii) = E_min*multiplicator**(iii - 1)
      end do

      max_error = 0.0d0

      ! loop over all tau time points
      do jquad = 1, num_integ_points

         current_point = tau_tj(jquad)

         ! y=exp(-x*|tau_k|)
         do iii = 1, num_x_nodes
            y_values(iii) = exp(-x_values(iii)*current_point)
         end do

         ! calculate mat_A
         do jjj = 1, num_integ_points
            do iii = 1, num_x_nodes
               omega = omega_tj(jjj)
               x_value = x_values(iii)
               mat_A(iii, jjj) = cos(current_point*omega)*2.0d0*x_value/(x_value**2 + omega**2)
            end do
         end do

         ! Singular value decomposition of mat_A
         call dgesdd('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
            _REGISTER_EXC("DGESDD returned ierr != 0")
            return
         end if

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         do jjj = 1, num_integ_points
            do iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            end do
         end do

         ! 2) U^T y
         call dgemm('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0d0, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0d0, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         call dgemm('N', 'N', num_integ_points, 1, num_x_nodes, 1.0d0, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0d0, weights_work, num_integ_points)

         weights(jquad, :) = weights_work(:)

         call calc_max_error_fit_omega_grid_with_cosine(max_error, current_point, omega_tj, weights_work, x_values, &
                                                        y_values, num_integ_points, num_x_nodes)

      end do ! jquad

      deallocate (x_values, y_values, mat_A, weights_work, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

   end subroutine get_l_sq_wghts_cos_tf_w_to_t


   pure subroutine calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                           y_values, num_integ_points, num_x_nodes)

      real(kind=dp), intent(inout)                       :: max_error, omega
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj, tau_wj_work, x_values, y_values
      integer, intent(in)                                :: num_integ_points, num_x_nodes

      integer                                            :: iii,kkk
      real(kind=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0d0

      do kkk = 1, num_x_nodes

         func_val = 0.0d0

         !call eval_fit_func_tau_grid_cosine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         ! calculate value of the fit function
         do iii = 1, num_integ_points
            func_val = func_val + tau_wj_work(iii)*cos(omega*tau_tj(iii))*exp(-x_values(kkk)*tau_tj(iii))
         end do

         
         if (abs(y_values(kkk) - func_val) > max_error_tmp) then
            max_error_tmp = abs(y_values(kkk) - func_val)
            func_val_temp = func_val
         end if

      end do

      if (max_error_tmp > max_error) then

         max_error = max_error_tmp

      end if

   end subroutine calc_max_error_fit_tau_grid_with_cosine


   pure subroutine eval_fit_func_tau_grid_cosine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

           !no necessary!
      real(kind=dp), intent(out)                         :: func_val
      real(kind=dp), intent(in)                          :: x_value
      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj, tau_wj_work
      real(kind=dp), intent(in)                          :: omega

      integer                                            :: iii

      func_val = 0.0d0

      do iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val + tau_wj_work(iii)*cos(omega*tau_tj(iii))*exp(-x_value*tau_tj(iii))

      end do

   end subroutine eval_fit_func_tau_grid_cosine


! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for sine transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   pure subroutine eval_fit_func_tau_grid_sine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      real(kind=dp), intent(inout)                       :: func_val
      real(kind=dp), intent(in)                          :: x_value
      integer, INTENT(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj, tau_wj_work
      real(kind=dp), intent(in)                          :: omega

      integer                                            :: iii

      func_val = 0.0d0

      do iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val + tau_wj_work(iii)*sin(omega*tau_tj(iii))*exp(-x_value*tau_tj(iii))

      end do

   end subroutine eval_fit_func_tau_grid_sine


   pure subroutine calc_max_error_fit_tau_grid_with_sine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                         y_values, num_integ_points, num_x_nodes)

      real(kind=dp), intent(inout)                       :: max_error, omega
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj, tau_wj_work, x_values, y_values
      integer, intent(in)                                :: num_integ_points, num_x_nodes

      integer                                            :: kkk
      real(kind=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0d0

      do kkk = 1, num_x_nodes

         func_val = 0.0d0

         call eval_fit_func_tau_grid_sine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         if (abs(y_values(kkk) - func_val) > max_error_tmp) then
            max_error_tmp = abs(y_values(kkk) - func_val)
            func_val_temp = func_val
         end if

      end do

      if (max_error_tmp > max_error) then

         max_error = max_error_tmp

      end if

   end subroutine calc_max_error_fit_tau_grid_with_sine


   pure subroutine eval_fit_func_omega_grid_cosine(func_val, x_value, num_integ_points, omega_tj, omega_wj_work, tau)
      real(kind=dp), intent(out)                         :: func_val
      real(kind=dp), intent(in)                          :: x_value
      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: omega_tj, omega_wj_work
      real(kind=dp), intent(in)                          :: tau

      integer                                            :: iii
      real(kind=dp)                                      :: omega

      func_val = 0.0d0

      do iii = 1, num_integ_points

         ! calculate value of the fit function
         omega = omega_tj(iii)
         func_val = func_val + omega_wj_work(iii)*cos(tau*omega)*2.0d0*x_value/(x_value**2 + omega**2)

      end do

   end subroutine eval_fit_func_omega_grid_cosine

   subroutine calc_max_error_fit_omega_grid_with_cosine(max_error, tau, omega_tj, omega_wj_work, x_values, &
                                                        y_values, num_integ_points, num_x_nodes)

      real(kind=dp), intent(inout)                       :: max_error
      real(kind=dp), intent(in)                          :: tau
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                         :: omega_tj, omega_wj_work, x_values, &
                                                               y_values
      integer, intent(in)                                :: num_integ_points, num_x_nodes

      integer                                            :: kkk
      real(kind=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0d0

      do kkk = 1, num_x_nodes

         func_val = 0.0d0

         call eval_fit_func_omega_grid_cosine(func_val, x_values(kkk), num_integ_points, omega_tj, omega_wj_work, tau)

         if (abs(y_values(kkk) - func_val) > max_error_tmp) then
            max_error_tmp = abs(y_values(kkk) - func_val)
            func_val_temp = func_val
         end if

      end do

      if (max_error_tmp > max_error) then

         max_error = max_error_tmp

      end if

   end subroutine calc_max_error_fit_omega_grid_with_cosine

   !> TODO(Maryam) Remove temporary routine
   subroutine get_minimax_grid(num_integ_points, &
           tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
           weights_cos_tf_w_to_t, weights_sin_tf_t_to_w,ounit)

      integer, parameter                                 :: unit_nr = 1
      integer, intent(in), optional                      :: ounit
      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
              intent(out)                                     :: tau_tj, tau_wj
      real(kind=dp), intent(in), optional                :: a_scaling
      real(kind=dp), intent(out), optional               :: e_fermi
      real(kind=dp), allocatable, dimension(:), &
              intent(inout), optional                         :: tj, wj
      real(kind=dp), allocatable, dimension(:, :), &
              intent(out), optional                           :: weights_cos_tf_t_to_w, &
              weights_cos_tf_w_to_t, &
              weights_sin_tf_t_to_w

      integer, parameter                                 :: num_points_per_magnitude = 200

      integer                                            :: ierr, jquad, i_exp
      real(kind=dp)                                      :: E_Range, Emax, Emin, max_error_min, &
              scaling, Range_from_i_exp

      real(kind=dp), allocatable, dimension(:)           :: x_tw

      ierr = 0
      allocate (tj(num_integ_points))
      allocate (wj(num_integ_points))
      allocate (tau_tj(0:num_integ_points))
      allocate (tau_wj(num_integ_points))
      allocate (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
      allocate (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
      allocate (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
      allocate (x_tw(2*num_integ_points))
      x_tw = 0.0d0
      ! @Maryam: Start of print statements that have been added

      do i_exp = 0, 50

         Range_from_i_exp = 1.58489319246d0**i_exp

         !IF (unit_nr > 0) then
         if (present(ounit)) then
            ! GET AND PRINT FREQUENCY GRIDS

            call get_rpa_minimax_grids(num_integ_points, Range_from_i_exp, x_tw, ierr)
            stop "The grid size you choose is not available."

            write (ounit, fmt="(T3,A,T66,F15.4)") "Range for the minimax approximation:", Range_from_i_exp
            write (ounit, fmt="(T3,A)") "minimax frequency omega_i     weight of frequency omega_i"
            do jquad = 1, num_integ_points
               write (ounit, fmt="(T15,F20.10,F20.10)") x_tw(jquad), x_tw(jquad + num_integ_points)
            end do

            ! GET AND PRINT TIME GRIDS

            if (num_integ_points .LE. 5) then
               write(*,*)"The grid size you choose is not available."
            else
               call get_exp_minimax_coeff_gw(num_integ_points, Range_from_i_exp, x_tw, ierr)
            end if

            write (ounit, fmt="(T3,A)") "minimax time tau_j     weight of time tau_j"
            do jquad = 1, num_integ_points
               write (ounit, fmt="(T15,F20.10,F20.10)") x_tw(jquad), x_tw(jquad + num_integ_points)
            end do

            write (ounit, fmt="(T3,A)") " "


         end if

      end do

      ! @Maryam: End of print statements that have been added

      deallocate (x_tw)

      do i_exp = 0, 50

         Range_from_i_exp = 1.58489319246d0**i_exp

         Emin = 1.0d0
         Emax = Range_from_i_exp
         E_Range = Emax/Emin



         !         if (.NOT. do_ri_sos_laplace_mp2) then
         allocate (x_tw(2*num_integ_points))
         x_tw = 0.0d0
         tj = 0.0d0
         wj = 0.0d0
         tau_tj = 0.0d0
         tau_wj = 0.0d0
         weights_cos_tf_t_to_w = 0.0d0
         weights_cos_tf_w_to_t = 0.0d0
         weights_sin_tf_t_to_w = 0.0d0

         call get_rpa_minimax_grids(num_integ_points, E_Range, x_tw, ierr)
         stop "The grid size you choose is not available."

         do jquad = 1, num_integ_points
            tj(jquad) = x_tw(jquad)
            wj(jquad) = x_tw(jquad + num_integ_points)
         end do

         deallocate (x_tw)



         ! scale the minimax parameters
         tj(:) = tj(:)*Emin
         wj(:) = wj(:)*Emin


         ! set up the minimax time grid
         !     if (do_im_time .OR. do_ri_sos_laplace_mp2) then
         allocate (x_tw(2*num_integ_points))
         x_tw = 0.0d0


         if (num_integ_points .LE. 5) then
            write(*,*)"The grid size you choose is not available."
         else
            call get_exp_minimax_coeff_gw(num_integ_points, E_Range, x_tw, ierr)
         end if

         ! For RPA we include already a factor of two (see later steps)
         scaling = 2.0d0
         !         if (do_ri_sos_laplace_mp2) scaling = 1.0d0



         do jquad = 1, num_integ_points
            tau_tj(jquad) = x_tw(jquad)/scaling
            tau_wj(jquad) = x_tw(jquad + num_integ_points)/scaling
         end do

         deallocate (x_tw)

         ! scale grid from [1,R] to [Emin,Emax]
         !tau_tj(:) = tau_tj(:)/Emin*Emin
         !tau_wj(:) = tau_wj(:)/Emin*Emin

         call get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, tj, weights_cos_tf_t_to_w, &
                 Emin, Emax, max_error_min, num_points_per_magnitude, ierr)

         call get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, tj, weights_cos_tf_w_to_t, &
                 Emin, Emax, max_error_min, num_points_per_magnitude, ierr)

         !write(1,*) "I am writing tj", tj, "***************************"


         call get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, tj, weights_sin_tf_t_to_w, &
                 Emin, Emax, max_error_min, num_points_per_magnitude, ierr)

         if (present(ounit)) then
            write(ounit, fmt="(T3,A,T66,F15.4)") "Range for the minimax approximation:", Range_from_i_exp
            write(ounit,*)Range_from_i_exp, "weights_cos_tf_t_to_w", weights_cos_tf_t_to_w
            write(ounit,*)Range_from_i_exp, "weights_cos_tf_w_to_t", weights_cos_tf_w_to_t
            write(ounit,*)Range_from_i_exp, "weights_sin_tf_t_to_w", weights_sin_tf_t_to_w
         end if
         !IF (unit_nr > 0) then
         !WRITE (UNIT=unit_nr, fmt="(T3,A,T66,ES15.2)") &
         ! "MINIMAX_INFO| Maximum deviation of the imag. time fit:", max_error_min
         !END if

      end do

      !deallocate(tj, wj, tau_tj, tau_wj, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)

   end subroutine get_minimax_grid

end module mp2_grids
