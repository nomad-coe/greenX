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
  use kinds,           only: dp
  use error_handling,  only: register_exc
  use constants,       only: pi
  use minimax_tau,     only: get_acoef_weight_tau
  use minimax_rpa,     only: get_rpa_minimax_grids

  implicit none

  private

  !> Main entry point for client code.
  public :: gx_minimax_grid

contains

  ! **************************************************************************************************
  ! \brief Compute minimax grid for RPA energy and GW calculation on imaginary time/frequency domain.
  ! o num_points: Number of mesh points.
  ! o e_min: Minimum transition energy. Arbitrary units as we only need e_max/e_min
  ! o e_max: Maximum transition energy.
  ! o tau_points: imaginary time grid points
  ! o tau_weights: weights for imaginary time grid points weights
  ! o omega_points: imaginary frequency grid points
  ! o omega_weights: weights for imaginary frequency grid points
  ! o cosft_wt: weights for tau -> omega cosine transform. cos(w*t) factor is included.
  ! o cosft_tw: weights for omega -> tau cosine transform. cos(w*t) factor is included.
  ! o sinft_wt: weights for tau -> omega sine transform. sin(w*t) factor is included.
  ! o max_errors: Max error for the three kind of transforms (same order as previous args)
  ! o cosft_duality_error. Max_{ij} |AB - I| where A and B are the cosft_wt and cosft_tw matrices.
  ! o ierr: Exit status
  ! **************************************************************************************************
  subroutine gx_minimax_grid(num_points, e_min, e_max, &
       tau_points, tau_weights, omega_points, omega_weights, &
       cosft_wt, cosft_tw, sinft_wt, &
       max_errors, cosft_duality_error, ierr)

    integer, intent(in)                               :: num_points
    real(dp), intent(in)                              :: e_min, e_max
    real(dp), allocatable, dimension(:), &  
         intent(out)                                  :: tau_points, tau_weights
    real(dp), allocatable, dimension(:), &
         intent(out)                                  :: omega_points(:), omega_weights(:)
    real(dp), allocatable, dimension(:, :), &
         intent(out)                                  :: cosft_wt(:, :), cosft_tw(:, :), sinft_wt(:, :)
    real(dp), intent(out)                             :: max_errors(3), cosft_duality_error
    integer, intent(out)                              :: ierr

    integer, parameter                                :: cos_t_to_cos_w = 1
    integer, parameter                                :: cos_w_to_cos_t = 2
    integer, parameter                                :: sin_t_to_sin_w = 3
    integer                                           :: i_point, j_point
    real(dp)                                          :: e_range, scaling
    real(dp), allocatable                             :: x_tw(:), mat(:, :)

    e_range = e_max/e_min
    ierr = 0

    allocate (x_tw(2*num_points))

    call get_rpa_minimax_grids(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    allocate (omega_points(num_points))
    allocate (omega_weights(num_points))

    do i_point = 1, num_points
       omega_points(i_point) = x_tw(i_point)
       omega_weights(i_point) = x_tw(i_point + num_points)
    end do

    ! scale the minimax parameters
    omega_points(:) = omega_points(:)*e_min
    omega_weights(:) = omega_weights(:)*e_min

    ! set up the minimax time grid
    call get_acoef_weight_tau(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    ! For RPA we include already a factor of two (see later steps)
    scaling = 2.0d0

    allocate (tau_points(num_points))
    allocate (tau_weights(num_points))

    do i_point = 1, num_points
       tau_points(i_point) = x_tw(i_point)/scaling
       tau_weights(i_point) = x_tw(i_point + num_points)/scaling
    end do

    ! scale grid from [1,R] to [e_min,e_max]
    tau_points(:) = tau_points(:)/e_min
    tau_weights(:) = tau_weights(:)/e_min

    allocate (cosft_wt(num_points, num_points))
    allocate (cosft_tw(num_points, num_points))
    allocate (sinft_wt(num_points, num_points))

    ! get the weights for the cosine transform W^c(it) -> W^c(iw)
    call get_transformation_weights(num_points, tau_points, omega_points, cosft_wt, e_min, e_max, &
         max_errors(1), cos_t_to_cos_w, ierr)
    if (ierr /= 0) return

    ! get the weights for the cosine transform W^c(iw) -> W^c(it)
    call get_transformation_weights(num_points, tau_points, omega_points, cosft_tw, e_min, e_max, &
         max_errors(2), cos_w_to_cos_t, ierr)
    if (ierr /= 0) return

    ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
    call get_transformation_weights(num_points, tau_points, omega_points, sinft_wt, e_min, e_max, &
         max_errors(3), sin_t_to_sin_w, ierr)
    if (ierr /= 0) return

    ! Compute the actual weights used for the inhomogeneous cosine/ FT and check whether
    ! the two matrices for the forward/backward transform are the inverse of each other.
    do i_point = 1, num_points
       do j_point = 1, num_points
          cosft_wt(j_point, i_point) = cosft_wt(j_point, i_point)*cos(tau_points(i_point)*omega_points(j_point))
          cosft_tw(i_point, j_point) = cosft_tw(i_point, j_point)*cos(tau_points(i_point)*omega_points(j_point))
          sinft_wt(j_point, i_point) = sinft_wt(j_point, i_point)*sin(tau_points(i_point)*omega_points(j_point))
       end do
    end do

    allocate (mat(num_points, num_points))
    mat(:, :) = matmul(cosft_wt, cosft_tw)
    do i_point = 1, num_points
       mat(i_point, i_point) = mat(i_point, i_point) - 1.0d0
    end do
    cosft_duality_error = maxval(abs(mat))

    deallocate (mat)
    deallocate (x_tw)

  end subroutine gx_minimax_grid

  ! **************************************************************************************************
  ! /brief get the weights eiter for the cosine/sin transformation for tau to omega or viceversa
  ! o num_points: Number of mesh points.
  ! o tau_points: imaginary time grid points
  ! o omega_points: imaginary frequency grid points
  ! o weights: corresponding tranformation weights
  ! o e_min: Minimum transition energy.
  ! o e_max: Maximum transition energy.
  ! o max_errors: Max error for the three kind of transforms
  ! o transformation type : 1 the cosine transform cos(it) -> cos(iw)
  !                       : 2 the cosine transform cos(iw) -> cos(it)
  !                       : 3 the sine transform   sin(it) -> sin(iw)
  ! o ierr: exit status
  ! **************************************************************************************************
  subroutine get_transformation_weights(num_points, tau_points, omega_points, weights, e_min, e_max, &
       max_error, transformation_type, ierr)

    integer, intent(in)                                :: num_points
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                    :: tau_points
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                    :: omega_points
    real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                 :: weights
    real(kind=dp), intent(in)                          :: e_min, e_max
    real(kind=dp), intent(inout)                       :: max_error
    integer, intent(in)                                :: transformation_type
    integer, intent(out)                               :: ierr

    integer                                            :: i_node, i_point, j_point, k_point, &
         num_x_nodes
    integer, parameter                                 :: nodes_factor = 200
    real(kind=dp)                                      :: current_point, x_factor
    real(kind=dp), allocatable, dimension(:)           :: weights_work, x_mu, psi
    real(kind=dp), allocatable, dimension(:, :)        :: mat_A

    integer                                            :: lwork
    integer, allocatable, dimension(:)                 :: iwork      
    real(kind=dp), allocatable, dimension(:)           :: vec_S, vec_UT_psi, work     
    real(kind=dp), allocatable, dimension(:, :)        :: mat_U, mat_VT, mat_VT_s

    ierr = 0

    allocate (weights_work(num_points))
    weights_work = 0.0d0

    ! compute the number of x nodes per magnitude points per 10-interval
    num_x_nodes = (int(log10(e_max/e_min)) + 1)*nodes_factor

    ! make sure that the number of x nodes are at least as many integration points
    num_x_nodes = max(num_x_nodes, num_points)

    allocate (x_mu(num_x_nodes))
    x_mu = 0.0d0
    allocate (psi(num_x_nodes))
    psi = 0.0d0

    ! Allocations for the BLAS routines
    ! double the value nessary for 'A' to achieve good performance
    lwork = 8*num_points*num_points + 12*num_points + 2*num_x_nodes
    allocate (iwork(8*num_points))
    iwork = 0
    allocate (work(lwork))
    work = 0.0d0      

    allocate (mat_A(num_x_nodes, num_points))
    mat_A = 0.0d0      
    allocate (mat_U(num_x_nodes, num_x_nodes))
    mat_U = 0.0d0
    allocate (mat_VT(num_x_nodes, num_points))
    mat_VT = 0.0d0
    allocate (mat_VT_s(num_points, num_x_nodes))
    mat_VT_s = 0.0d0
    allocate (vec_S(num_points))
    vec_S = 0.0d0      
    allocate (vec_UT_psi(num_x_nodes))
    vec_UT_psi = 0.0d0

    ! set the x-mu logarithmically in the interval [e_min,e_max]
    x_factor = (e_max/e_min)**(1.0d0/(real(num_x_nodes, kind=dp) - 1.0d0))
    do i_node = 1, num_x_nodes
       x_mu(i_node) = e_min*x_factor**(i_node - 1)
    end do

    current_point = 0.0d0
    max_error = 0.0d0

    ! loop over all grid points
    do i_point = 1, num_points
       ! calculate psi and mat_A
       call calculate_psi_and_mat_A(num_points, tau_points, omega_points, num_x_nodes, x_mu, psi, &
            mat_A, i_point, current_point, transformation_type)

       ! Singular value decomposition of mat_A = U*Sigma*V^T
       call dgesdd('A', num_x_nodes, num_points, mat_A, num_x_nodes, vec_S, mat_U, num_x_nodes, &
            mat_VT, num_x_nodes, work, lwork, iwork, ierr)

       if (ierr /= 0) then
          _REGISTER_EXC("DGESDD returned ierr != 0")
          return
       end if

       ! integration weights = (V Sigma^-1 U^T)*psi

       ! 1) V * Sigma^-1
       do j_point = 1, num_points
          do k_point = 1, num_points
             mat_VT_s(k_point, j_point) = mat_VT(j_point, k_point)/vec_S(j_point)
          end do ! k_point
       end do ! j_point

       ! 2) (U^T)*psi
       call dgemm('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0d0, mat_U, num_x_nodes, psi, num_x_nodes, &
            0.0d0, vec_UT_psi, num_x_nodes)

       ! 3) (V*Sigma^-1) * (U^T*psi)
       call dgemm('N', 'N', num_points, 1, num_x_nodes, 1.0d0, mat_VT_s, num_points, vec_UT_psi, &
            num_x_nodes, 0.0d0, weights_work, num_points)

       weights(i_point, :) = weights_work(:)

       ! calculate the maximum error of the fitting
       call calculate_max_error(num_points, tau_points, omega_points, weights_work, num_x_nodes, x_mu, &
            psi, current_point, max_error, transformation_type)
    end do ! i_point

    deallocate (x_mu, psi, mat_A, weights_work, vec_S, mat_U, mat_VT, work, iwork, mat_VT_s, vec_UT_psi)

  end subroutine get_transformation_weights

  ! **************************************************************************************************
  ! /brief calculate the auxiliary matrix for cosine/sin transformation for tau to omega or viceversa
  ! o num_points: Number of mesh points.
  ! o tau_points: imaginary time grid points
  ! o omega_points: imaginary frequency grid points
  ! o num_x_nodes: Number of node in the interval [e_min,e_max]
  ! o x_mu : Transition energy (nodes in the interval [e_min,e_max])
  ! o psi: corresponding auxiliary function (see transformation type definition)
  ! o mat_A: auxiliary matrix (see transformation type definition)  
  ! o i_point: pointer for the current grid point
  ! o current_point:  current grid point ether omega(i_point) or tau_(i_point)
  ! o transformation type :
  !   (1) the cosine transform cos(it) -> cos(iw): psi(omega,x), mat_A = cos(omega*tau)*psi(tau,x) 
  !   (2) the cosine transform cos(iw) -> cos(it): psi(tau,x)  , mat_A = cos(omega*tau)*psi(omega,x)
  !   (3) the sine transform   sin(it) -> sin(iw): psi(omega,x), mat_A = sin(omega*tau)*psi(tau,x)
  ! **************************************************************************************************
  subroutine calculate_psi_and_mat_A(num_points, tau_points, omega_points, num_x_nodes, x_mu, psi, &
       mat_A, i_point, current_point, transformation_type)

    integer, intent(in)                                :: num_points, num_x_nodes, i_point
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                    :: tau_points, omega_points, x_mu
    real(kind=dp), allocatable, dimension(:), &          
         intent(inout)                                 :: psi         
    real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                 :: mat_A
    real(kind=dp), intent(inout)                       :: current_point
    integer, intent(in)                                :: transformation_type

    integer                                            :: i_node, j_point
    real(kind=dp)                                      :: tau, omega

    ! the cosine transform cos(it) -> cos(iw)
    if (transformation_type.eq.1) then
       omega = omega_points(i_point)
       current_point = omega

       ! psi(omega_k,x) = 2x/(x^2+omega_k^2)
       do i_node = 1, num_x_nodes
          psi(i_node) = 2.0d0*x_mu(i_node)/((x_mu(i_node))**2 + omega**2)
       end do

       ! mat_A = cos(omega_k * tau) psi(tau,x)
       do j_point = 1, num_points
          tau = tau_points(j_point)
          do i_node = 1, num_x_nodes
             mat_A(i_node, j_point) = cos(omega*tau)*exp(-x_mu(i_node)*tau)
          end do
       end do

    ! the cosine transform cos(iw) -> cos(it)  
    else if (transformation_type.eq.2) then
       tau = tau_points(i_point)
       current_point = tau

       ! psi(tau_k,x) = =exp(-x*|tau_k|)
       do i_node = 1, num_x_nodes
          psi(i_node) = exp(-x_mu(i_node)*tau)
       end do

       ! mat_A = cos(tau_k,omega) psi(omega,x)
       do j_point = 1, num_points
          omega = omega_points(j_point)
          do i_node = 1, num_x_nodes
             mat_A(i_node, j_point) = cos(tau*omega)*2.0d0*x_mu(i_node)/(x_mu(i_node)**2 + omega**2)
          end do
       end do

    ! the sine transform sin(it) -> sin(iw)         
    else if (transformation_type.eq.3) then
       omega = omega_points(i_point)
       current_point = omega

       ! psi(omega_k,x) = 2*omega_k/(x^2+omega_k^2)
       do i_node = 1, num_x_nodes
          psi(i_node) = 2.0d0*omega/((x_mu(i_node))**2 + omega**2)
       end do

       ! mat_A = sin(omega_k,tau)*phi(tau,x)
       do j_point = 1, num_points
          tau = tau_points(j_point)
          do i_node = 1, num_x_nodes
             mat_A(i_node, j_point) = sin(omega*tau)*exp(-x_mu(i_node)*tau)
          end do
       end do
    end if

  end subroutine calculate_psi_and_mat_A

  ! **************************************************************************************************
  ! /brief calculate the error of the fit function for the cosine/sin transformation for tau to omega or viceversa
  ! o num_points: Number of mesh points.
  ! o tau_points: imaginary time grid points
  ! o omega_points: imaginary frequency grid points
  ! o weights_work: work vector for the transformation weights
  ! o num_x_nodes: Number of node in the interval [e_min,e_max]
  ! o x_mu : Transition energy (nodes in the interval [e_min,e_max])
  ! o psi: corresponding auxiliary function (see transformation type definition)
  ! o current_point:  current grid point ether omega(i_point) or tau_(i_point)
  ! o max_errors: Max error for the three kind of transforms   
  ! o transformation type : 1 fit function for the cosine transform cos(it) -> cos(iw), psi(omeaga,x)
  !                       : 2 fit function for the cosine transform cos(iw) -> cos(it), psi(tau,x)
  !                       : 3 fit function for the sine transform   sin(it) -> sin(iw), psi(omega,x)
  ! **************************************************************************************************
  subroutine calculate_max_error(num_points, tau_points, omega_points, weights_work, num_x_nodes, x_mu, &
       psi, current_point, max_error, transformation_type)

    real(kind=dp), intent(out)                       :: max_error
    real(kind=dp), intent(in)                        :: current_point      
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                  :: tau_points, omega_points, x_mu, psi, &
         weights_work
    integer, intent(in)                              :: num_points, num_x_nodes
    integer, intent(in)                              :: transformation_type      

    integer                                          :: i_node,i_point
    real(kind=dp)                                    :: func_val, func_val_temp, max_error_tmp, &
         tau, omega, x_val

    max_error_tmp = 0.0d0

    ! the cosine transform cos(it) -> cos(iw)
    if (transformation_type.eq.1) then
       omega=current_point

       do i_node = 1, num_x_nodes
          func_val = 0.0d0
          ! calculate value of the fit function f(x) = f(x) - weights(omega)cos(omega*tau)psi(tau.x)
          do i_point = 1, num_points
             tau = tau_points(i_point) 
             func_val = func_val + weights_work(i_point)*cos(omega*tau)*exp(-x_mu(i_node)*tau)
          end do

          if (abs(psi(i_node) - func_val) > max_error_tmp) then
             max_error_tmp = abs(psi(i_node) - func_val)
             func_val_temp = func_val
          end if
       end do

    ! the cosine transform cos(iw) -> cos(it)
    else if (transformation_type.eq.2) then
       tau = current_point

       do i_node = 1, num_x_nodes
          func_val = 0.0d0
          x_val=x_mu(i_node)
          ! calculate value of the fit function f(x) = f(x) - weights(tau)cos(omega*tau)psi(omega.x)
          do i_point = 1, num_points
             omega = omega_points(i_point)
             func_val = func_val +  weights_work(i_point)*cos(tau*omega)*2.0d0*x_val/(x_val**2 + omega**2)
          end do

          if (abs(psi(i_node) - func_val) > max_error_tmp) THEN
             max_error_tmp = abs(psi(i_node) - func_val)
             func_val_temp = func_val
          end if
       end do

    ! the sine transform sin(it) -> sin(iw)
    else if (transformation_type.eq.3) then
       omega = current_point

       do i_node = 1, num_x_nodes
          func_val = 0.0d0
          ! calculate value of the fit function f(x) = f(x) - weights(omega)sin(omega*tau)psi(tau.x)
          do i_point = 1, num_points
             tau = tau_points(i_point)
             func_val = func_val +  weights_work(i_point)*sin(omega*tau)*exp(-x_mu(i_node)*tau)
          end do

          if (abs(psi(i_node) - func_val) > max_error_tmp) then
             max_error_tmp = abs(psi(i_node) - func_val)
             func_val_temp = func_val
          end if
       end do
    end if

    if (max_error_tmp > max_error) then
       max_error = max_error_tmp
    end if

  end subroutine calculate_max_error

end module minimax_grids
