!  Copyright (C) 2020-2022 Green-X library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief Routines to calculate frequency and time grids (integration points and weights)
!> for correlation methods as well as weights for the inhomogeneous cosine/sine transform.
!>
!> NB: When dealing with runtime exceptions, we set ierr to a non-zero value and return immediately
!  to the caller so we don't need to goto to a cleanup section at the end of the procedure.
!  Assume -std=f2008: i.e. allocatable arrays are automatically deallocated when going out of scope.
! **************************************************************************************************

module minimax_grids
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
!> \param[out] tau_weights: weights for tau mesh (tau_wj in CP2K)
!> \param[out] omega_mesh: imaginary frequency mesh (tj in CP2K)
!> \param[out] omega_weights: weights for imaginary frequency mesh (wj in CP2K)
!> \param[out] cosft_wt: weights for tau -> w cosine transform. cos(w*t) factor is included.
!> \param[out] cosft_tw: weights for w -> tau cosine transform. cos(w*t) factor is included.
!> \param[out] sinft_wt: weights for tau -> w sine transform. sin(w*t) factor is included.
!> \param[out] max_errors: Max error for the three kind of transforms (same order as previous args)
!> \param[out] cosft_duality_error. Max_{ij} |AB - I| where A and B are the cosft_wt and cosft_tw matrices.
!> \param[out] ierr: Exit status
   subroutine gx_minimax_grid(num_integ_points, emin, emax, &
                              tau_mesh, tau_weights, omega_mesh, omega_weights, &
                              cosft_wt, cosft_tw, sinft_wt, &
                              max_errors, cosft_duality_error, ierr)

      integer, intent(in) :: num_integ_points
      real(dp), intent(in) :: emin, emax
      real(dp), allocatable, intent(out) :: tau_mesh(:), tau_weights(:)
      real(dp), allocatable, intent(out) :: omega_mesh(:), omega_weights(:)
      real(dp), allocatable, intent(out) :: cosft_wt(:, :), cosft_tw(:, :), sinft_wt(:, :)
      real(dp), intent(out) :: max_errors(3), cosft_duality_error
      integer, intent(out) :: ierr

      integer, parameter                                 :: num_points_per_magnitude = 200
      integer, parameter                                 :: cos_t_to_cos_w = 1
      integer, parameter                                 :: cos_w_to_cos_t = 2
      integer, parameter                                 :: sin_t_to_sin_w = 3
      integer                                            :: i_point, j_point
      real(dp)                                           :: e_range, scaling
      real(dp), allocatable                              :: x_tw(:), mat(:, :)

      E_Range = Emax/Emin
      ierr = 0

      allocate (x_tw(2*num_integ_points))

      call get_rpa_minimax_grids(num_integ_points, E_Range, x_tw, ierr)
      if (ierr /= 0) return

      allocate (omega_mesh(num_integ_points))
      allocate (omega_weights(num_integ_points))

      do i_point = 1, num_integ_points
         omega_mesh(i_point) = x_tw(i_point)
         omega_weights(i_point) = x_tw(i_point + num_integ_points)
      end do

      ! scale the minimax parameters
      omega_mesh(:) = omega_mesh(:)*Emin
      omega_weights(:) = omega_weights(:)*Emin

      ! set up the minimax time grid
      call get_exp_minimax_coeff_gw(num_integ_points, E_Range, x_tw, ierr)
      if (ierr /= 0) return

      ! For RPA we include already a factor of two (see later steps)
      scaling = 2.0d0

      allocate (tau_mesh(num_integ_points))
      allocate (tau_weights(num_integ_points))

      do i_point = 1, num_integ_points
         tau_mesh(i_point) = x_tw(i_point)/scaling
         tau_weights(i_point) = x_tw(i_point + num_integ_points)/scaling
      end do

      ! scale grid from [1,R] to [Emin,Emax]
      tau_mesh(:) = tau_mesh(:)/Emin
      tau_weights(:) = tau_weights(:)/Emin

      allocate (cosft_wt(num_integ_points, num_integ_points))
      allocate (cosft_tw(num_integ_points, num_integ_points))
      allocate (sinft_wt(num_integ_points, num_integ_points))

      ! get the weights for the cosine transform W^c(it) -> W^c(iw)
      call get_tranformation_weights(num_integ_points, tau_mesh, omega_mesh, cosft_wt, Emin, Emax, &
                                     max_errors(1), num_points_per_magnitude, cos_t_to_cos_w, ierr)
      if (ierr /= 0) return

      ! get the weights for the cosine transform W^c(iw) -> W^c(it)
      call get_tranformation_weights(num_integ_points, tau_mesh, omega_mesh, cosft_tw, Emin, Emax, &
                                     max_errors(2), num_points_per_magnitude, cos_w_to_cos_t, ierr)
      if (ierr /= 0) return

      ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
      call get_tranformation_weights(num_integ_points, tau_mesh, omega_mesh, sinft_wt, Emin, Emax, &
                                     max_errors(3), num_points_per_magnitude, sin_t_to_sin_w, ierr)
      if (ierr /= 0) return

      ! Compute the actual weights used for the inhomogeneous cosine/ FT and check whether
      ! the two matrices for the forward/backward transform are the inverse of each other.
      do i_point = 1, num_integ_points
         do j_point = 1, num_integ_points
            cosft_wt(j_point, i_point) = cosft_wt(j_point, i_point)*cos(tau_mesh(i_point)*omega_mesh(j_point))
            cosft_tw(i_point, j_point) = cosft_tw(i_point, j_point)*cos(tau_mesh(i_point)*omega_mesh(j_point))
            sinft_wt(j_point, i_point) = sinft_wt(j_point, i_point)*sin(tau_mesh(i_point)*omega_mesh(j_point))
         end do
      end do

      allocate (mat(num_integ_points, num_integ_points))
      mat(:, :) = matmul(cosft_wt, cosft_tw)
      do i_point = 1, num_integ_points
         mat(i_point, i_point) = mat(i_point, i_point) - 1.0d0
      end do
      cosft_duality_error = maxval(abs(mat))

      deallocate (mat)
      deallocate (x_tw)

   end subroutine gx_minimax_grid


   subroutine get_tranformation_weights(num_integ_points, tau_tj, omega_tj, weights,E_min, E_max, &
                                           max_error, num_points_per_magnitude, transformation_type, ierr)

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
      integer, intent(in)                                :: transformation_type
      integer, intent(out)                               :: ierr

      integer                                            :: i_node, i_point, j_point, k_point, &
                                                            lwork, num_x_nodes
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
      do i_node = 1, num_x_nodes
         x_values(i_node) = E_min*multiplicator**(i_node - 1)
      end do

      current_point = 0.0d0
      max_error = 0.0d0

      ! loop over all omega frequency points
      do i_point = 1, num_integ_points

        ! calculate mat_A
         call calculate_mat_A(num_integ_points, tau_tj, omega_tj, x_values, y_values, &
                              num_x_nodes, mat_A, i_point, current_point, transformation_type)

         ! Singular value decomposition of mat_A
         call dgesdd('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, ierr)

         if (ierr /= 0) then
            _REGISTER_EXC("DGESDD returned ierr != 0")
            return
         end if

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         do j_point = 1, num_integ_points
            do k_point = 1, num_integ_points
               mat_SinvVSinvSigma(k_point, j_point) = mat_SinvVSinvT(j_point, k_point)/sing_values(j_point)
            end do
         end do

         ! 2) U^T y
         call dgemm('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0d0, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0d0, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         call dgemm('N', 'N', num_integ_points, 1, num_x_nodes, 1.0d0, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0d0, weights_work, num_integ_points)

         weights(i_point, :) = weights_work(:)

         ! calculate the maximum error of the fitting
         call calculate_max_error(max_error, current_point, tau_tj, omega_tj, weights_work, x_values, &
                                  y_values, num_integ_points, num_x_nodes, transformation_type)                                      
      end do ! i_point

      deallocate (x_values, y_values, mat_A, weights_work, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

   end subroutine get_tranformation_weights

   subroutine calculate_mat_A(num_integ_points, tau_tj, omega_tj, x_values, y_values, &
                              num_x_nodes, mat_A, i_point, current_point, transformation_type)
 
      integer, intent(in)                                :: num_integ_points, num_x_nodes, i_point
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj, omega_tj, x_values
      real(kind=dp), allocatable, dimension(:), &          
         intent(inout)                                   :: y_values         
      real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                   :: mat_A
      real(kind=dp), intent(inout)                       :: current_point
      integer, intent(in)                                :: transformation_type

      integer                                            :: i_node, j_point
      real(kind=dp)                                      :: tau, omega


      if (transformation_type.eq.1) then

         omega = omega_tj(i_point)
         current_point = omega

         ! y=2x/(x^2+omega_k^2)
         do i_node = 1, num_x_nodes
            y_values(i_node) = 2.0d0*x_values(i_node)/((x_values(i_node))**2 + omega**2)
         end do

         ! calculate mat_A
         do j_point = 1, num_integ_points
            do i_node = 1, num_x_nodes
               mat_A(i_node, j_point) = cos(omega*tau_tj(j_point))*exp(-x_values(i_node)*tau_tj(j_point))
            end do
        end do

      else if (transformation_type.eq.2) then

         tau = tau_tj(i_point)
         current_point = tau

         ! y=exp(-x*|tau_k|)
         do i_node = 1, num_x_nodes
            y_values(i_node) = exp(-x_values(i_node)*tau)
         end do

         ! calculate mat_A
         do j_point = 1, num_integ_points
            omega = omega_tj(j_point)
            do i_node = 1, num_x_nodes
               mat_A(i_node, j_point) = cos(tau*omega)*2.0d0*x_values(i_node)/(x_values(i_node)**2 + omega**2)
            end do
         end do

      else if (transformation_type.eq.3) then

         omega = omega_tj(i_point)
         current_point = omega

         ! y=2*omega_k/(x^2+omega_k^2)
         do i_node = 1, num_x_nodes
            y_values(i_node) = 2.0d0*omega/((x_values(i_node))**2 + omega**2)
         end do

         ! calculate mat_A
         do j_point = 1, num_integ_points
            do i_node = 1, num_x_nodes
               mat_A(i_node, j_point) = sin(omega*tau_tj(j_point))*exp(-x_values(i_node)*tau_tj(j_point))
            end do
         end do

      end if

   end subroutine calculate_mat_A

   subroutine calculate_max_error(max_error, current_point, tau_tj, omega_tj, weights_work, x_values, &
                                  y_values, num_integ_points, num_x_nodes, transformation_type)

      real(kind=dp), intent(inout)                       :: max_error, current_point
      real(kind=dp), allocatable, dimension(:), &
         intent(in)                                      :: tau_tj, omega_tj, x_values, y_values, &
                                                            weights_work
      integer, intent(in)                                :: num_integ_points, num_x_nodes
      integer, intent(in)                                :: transformation_type      

      integer                                            :: i_node,i_point
      real(kind=dp)                                      :: func_val, func_val_temp, max_error_tmp, &
                                                            tau, omega, x_value

      max_error_tmp = 0.0d0

      !!!! Cosine t to w!!!
      if (transformation_type.eq.1) then

         omega=current_point

         do i_node = 1, num_x_nodes

            func_val = 0.0d0
            x_value=x_values(i_node)

            ! calculate value of the fit function
            do i_point = 1, num_integ_points
              tau = tau_tj(i_point) 
              func_val = func_val + weights_work(i_point)*cos(omega*tau)*exp(-x_value*tau)
            end do

            if (abs(y_values(i_node) - func_val) > max_error_tmp) then
               max_error_tmp = abs(y_values(i_node) - func_val)
               func_val_temp = func_val
            end if

         end do

      !!!! Cosine w to t  !!!!
     else if (transformation_type.eq.2) then

        tau = current_point

        do i_node = 1, num_x_nodes

           func_val = 0.0d0
           x_value=x_values(i_node)

           ! calculate value of the fit function
           do i_point = 1, num_integ_points
              omega = omega_tj(i_point)
              func_val = func_val +  weights_work(i_point)*cos(tau*omega)*2.0d0*x_value/(x_value**2 + omega**2)
           end do

           if (abs(y_values(i_node) - func_val) > max_error_tmp) THEN
              max_error_tmp = abs(y_values(i_node) - func_val)
              func_val_temp = func_val
           end if

        end do

      !!!! sine t to w !!!!
      else if (transformation_type.eq.3) then

         omega = current_point

         do i_node = 1, num_x_nodes

            func_val = 0.0d0
            x_value=x_values(i_node)

            ! calculate value of the fit function
            do i_point = 1, num_integ_points
               tau = tau_tj(i_point)
               func_val = func_val +  weights_work(i_point)*sin(omega*tau)*exp(-x_value*tau)
            end do

            if (abs(y_values(i_node) - func_val) > max_error_tmp) then
               max_error_tmp = abs(y_values(i_node) - func_val)
              func_val_temp = func_val
            end if
         end do   
      end if 

      if (max_error_tmp > max_error) then

         max_error = max_error_tmp

      end if   

   end subroutine calculate_max_error

   !> TODO(Maryam) Remove temporary routine
   subroutine get_minimax_grid(num_integ_points, &
           tau_tj, tau_wj, a_scaling, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
           weights_cos_tf_w_to_t, weights_sin_tf_t_to_w,ounit)

      integer, parameter                                 :: unit_nr = 1
      integer, intent(in), optional                      :: ounit
      integer, intent(in)                                :: num_integ_points
      real(kind=dp), allocatable, dimension(:), &
              intent(out)                                :: tau_tj, tau_wj
      real(kind=dp), intent(in), optional                :: a_scaling
      real(kind=dp), intent(out), optional               :: e_fermi
      real(kind=dp), allocatable, dimension(:), &
              intent(inout), optional                    :: tj, wj
      real(kind=dp), allocatable, dimension(:, :), &
              intent(out), optional                      :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w

      integer, parameter                                 :: num_points_per_magnitude = 200
      integer, parameter                                 :: cos_t_to_cos_w = 1
      integer, parameter                                 :: cos_w_to_cos_t = 2
      integer, parameter                                 :: sin_t_to_sin_w = 3      

      integer                                            :: ierr, i_point, i_exp
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
            do i_point = 1, num_integ_points
               write (ounit, fmt="(T15,F20.10,F20.10)") x_tw(i_point), x_tw(i_point + num_integ_points)
            end do

            ! GET AND PRINT TIME GRIDS

            if (num_integ_points .LE. 5) then
               write(*,*)"The grid size you choose is not available."
            else
               call get_exp_minimax_coeff_gw(num_integ_points, Range_from_i_exp, x_tw, ierr)
            end if

            write (ounit, fmt="(T3,A)") "minimax time tau_j     weight of time tau_j"
            do i_point = 1, num_integ_points
               write (ounit, fmt="(T15,F20.10,F20.10)") x_tw(i_point), x_tw(i_point + num_integ_points)
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

         do i_point = 1, num_integ_points
            tj(i_point) = x_tw(i_point)
            wj(i_point) = x_tw(i_point + num_integ_points)
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



         do i_point = 1, num_integ_points
            tau_tj(i_point) = x_tw(i_point)/scaling
            tau_wj(i_point) = x_tw(i_point + num_integ_points)/scaling
         end do

         deallocate (x_tw)

         ! scale grid from [1,R] to [Emin,Emax]
         !tau_tj(:) = tau_tj(:)/Emin*Emin
         !tau_wj(:) = tau_wj(:)/Emin*Emin

         call get_tranformation_weights(num_integ_points, tau_tj, tj, weights_cos_tf_t_to_w, Emin, Emax, &
                                        max_error_min, num_points_per_magnitude, cos_t_to_cos_w, ierr)

         call get_tranformation_weights(num_integ_points, tau_tj, tj, weights_cos_tf_w_to_t, Emin, Emax, &
                                        max_error_min, num_points_per_magnitude, cos_w_to_cos_t, ierr)

         call get_tranformation_weights(num_integ_points, tau_tj, tj, weights_sin_tf_t_to_w, Emin, Emax, &
                                        max_error_min, num_points_per_magnitude, sin_t_to_sin_w, ierr)

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

end module minimax_grids
