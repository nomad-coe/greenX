! **************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief This module contains the subroutines for the localized basis set component of the library
! ***************************************************************************************************
module localized_basis_environments

   use kinds,                 only: dp
   use error_handling,    only: register_exc   
   use localized_basis_types, only: separable_ri_types, &
                                    screened_coulomb_types
   use lapack_interfaces,     only: dgemm, dsyevx
   use gx_minimax,            only: gx_minimax_grid
   


   implicit none

   !*******************************************************************!  

   contains

  !> \brief Initialize the ri_rs types.
  !! @param[inout] ri_rs: type for the separable ri
   subroutine initialization(ri_rs)

   type(separable_ri_types) :: ri_rs

   ri_rs%basis%n_basbas=303
   ri_rs%basis%n_basis = 58 
   ri_rs%basis%n_loc_basbas=303
   ri_rs%basis%n_basis_pairs=1711
   ri_rs%n_points = 638

   ri_rs%error = 0.d0

   if (.not.allocated(ri_rs%ovlp2fn)) then
      allocate(ri_rs%ovlp2fn(ri_rs%basis%n_basis_pairs,ri_rs%n_points))
   end if

   if (.not.allocated(ri_rs%ovlp3fn)) then
      allocate(ri_rs%ovlp3fn(ri_rs%basis%n_basis_pairs,ri_rs%basis%n_basbas))
   end if

   if (.not.allocated(ri_rs%z_coeff)) then
      allocate(ri_rs%z_coeff(ri_rs%basis%n_basbas,ri_rs%n_points))
   end if
   
   end subroutine initialization

! *********************************************************************
!> brief allocate the minimax time-frequency grids 
!  o ri_rs --  separable resolution-of-the identity environment
! *********************************************************************
   subroutine initialize_minimax_grids(w_pq)

     type(screened_coulomb_types)               :: w_pq

     ! Hard coded (change later) 
     w_pq%minimax%n_points = 6

     if (.not.allocated(w_pq%minimax%cos_tf)) then
        allocate(w_pq%minimax%cos_tf(w_pq%minimax%n_points, w_pq%minimax%n_points))
     end if

     if (.not.allocated(w_pq%minimax%omega)) then
        allocate(w_pq%minimax%omega(w_pq%minimax%n_points))
     end if

     if (.not.allocated(w_pq%minimax%tau)) then
        allocate(w_pq%minimax%tau(w_pq%minimax%n_points))
     end if

     if (.not.allocated(w_pq%minimax%weights)) then
        allocate(w_pq%minimax%weights(w_pq%minimax%n_points))
     end if

   end subroutine initialize_minimax_grids
     

  !> \brief Deallocate the ri_rs types.
  !! @param[inout] ri_rs: type for the separable ri   
   subroutine deallocations(ri_rs)

   type(separable_ri_types) :: ri_rs

   if (allocated(ri_rs%ovlp2fn)) deallocate(ri_rs%ovlp2fn)

   if (allocated(ri_rs%ovlp3fn)) deallocate(ri_rs%ovlp3fn)

   if (allocated(ri_rs%z_coeff)) deallocate(ri_rs%z_coeff)   

   end subroutine deallocations

  !> \brief Compute the error between the RI-V and RI-RS three center overlap coefficients
  !! @param[in] ri_rs: Type for the separable ri
  !! @param[in] n_basis_pairs: Number of orbital basis pairs (dimension 1 of ovlp3fn array)
  !! @param[in] n_basbas: Number of auxiliary basis fuctions (dimension 2 of ovlp3fn array)
  !! @param[in] ovlp3fn: RI-V three-center overlap coefficients
  !! @param[out] n_basis_pairs: maximun error
   subroutine calculate_error(ri_rs, n_basis_pairs, n_basbas, ovlp3fn, error)

   integer n_basbas, n_basis_pairs
   real(kind=8)                                       :: error
   real(kind=8), dimension(n_basis_pairs,n_basbas)    :: ovlp3fn

   type(separable_ri_types) :: ri_rs   

   ! Local variables

   integer      :: i_basbas, i_pair
   real(kind=dp) :: tmp_error

   !*******************************************************************!

   do i_pair=1,n_basis_pairs
      do i_basbas=1, n_basbas
         tmp_error = abs(ri_rs%ovlp3fn(i_pair,i_basbas) - ovlp3fn(i_pair,i_basbas))
         error=max(error,tmp_error) 
      end do
   end do

   end subroutine calculate_error

  !> \brief Get the machine precision for safe threshold during the lapack diagonalization
  !! @param[in] safe_minimum: Type for the separable ri
   subroutine get_machine_precision(safe_minimum)

   real(kind=dp) :: safe_minimum
   real(kind=dp), external ::  dlamch

   safe_minimum = dlamch ('S') 

   end subroutine get_machine_precision

! **********************************************************************
!> brief get frequency-time minimax grids 
!  o w_pq --  Screened Coulomb environment
! **********************************************************************
  subroutine get_minimax_grids(w_pq)

    type(screened_coulomb_types)               :: w_pq

    !Local variables

    integer                                   :: info, n_points
    real(kind=8)                              :: e_tran_max, e_tran_min, duality_error
    real(kind=8), dimension(3)                :: max_errors
    real(kind=8), dimension(:),   allocatable :: freq_mesh, freq_weights
    real(kind=8), dimension(:),   allocatable :: tau_mesh, tau_weights
    real(kind=8), dimension(:,:), allocatable :: cos_tau_to_freq_weights
    real(kind=8), dimension(:,:), allocatable :: cos_freq_to_tau_weights
    real(kind=8), dimension(:,:), allocatable :: sinft_tau_to_freq_weights

    ! Initialization

    call initialize_minimax_grids(w_pq)    

    n_points = w_pq%minimax%n_points

    ! Get the transition energies

    ! call get_minimal_maximal_transition_energy(e_tran_min, e_tran_max)
      ! Hard coded for now
      e_tran_max = 0.260
      e_tran_min = 31.725

    if (e_tran_min .le. 0.d0) then
       call register_exc("Detected metal system, not supported in minimax grid")
       return
    end if

    ! Get imaginary time and frequency points and transformations matrices

    call gx_minimax_grid(n_points, e_tran_min, e_tran_max, &
         tau_mesh, tau_weights, freq_mesh, freq_weights, &
         cos_tau_to_freq_weights, cos_freq_to_tau_weights, &
         sinft_tau_to_freq_weights, max_errors, duality_error, info)
    if (info /= 0) then
        call register_exc("Error in getting the minimax grid")
       return
    end if

    w_pq%minimax%cos_tf(1: n_points,1: n_points) = cos_tau_to_freq_weights(1: n_points,1: n_points)
    w_pq%minimax%omega(1: n_points)= freq_mesh(1: n_points)
    w_pq%minimax%tau(1: n_points) = tau_mesh (1: n_points)
    w_pq%minimax%weights(1: n_points) = freq_weights(1: n_points)

  end subroutine get_minimax_grids
   

  !> \brief Compute the power of a matrix using the lapack diagonalization
  !! @param[inout] matrix: working array
  !! @param[in] n_dim: dimension of the working array
  !! @param[in] power: exponent of the power
  !! @param[in] threshold: threshold for the eigenvalue 
   subroutine power_genmat(matrix, n_dim, power, threshold)

   integer :: n_dim

   real(kind=dp)                              :: power, threshold
   real(kind=dp), dimension(n_dim,n_dim)      :: matrix

   ! Local variables

   integer                                   :: i_element, j_element, &
                                                n_nonsingular

   real(kind=dp)                              :: ev_sqrt
   real(kind=dp), parameter                   :: my_threshold = -1.d10
   real(kind=dp), dimension(:), allocatable   :: eigenvalues
   real(kind=dp), dimension(:,:), allocatable :: work


   allocate(eigenvalues(n_dim),work(n_dim,n_dim))
   eigenvalues=0._dp
   work=0._dp

   matrix = - matrix

   call diagonalize_genmat(matrix,n_dim,my_threshold,n_nonsingular,eigenvalues,work)

   eigenvalues = - eigenvalues   

   if (n_nonsingular /= n_dim) then
       write(*,*) n_dim, n_nonsingular
       write(*,*) eigenvalues
       write(*,*) "Unphysical eigenvalues"
       stop
   end if

   do i_element = 1, n_dim
      if (eigenvalues(i_element) <= threshold) exit
   end do
   n_nonsingular = min(i_element-1, n_dim)

   do i_element = 1, n_nonsingular
      ev_sqrt = sqrt(eigenvalues(i_element))
      do j_element = 1, n_dim, 1
         work(j_element, i_element) = &
         work(j_element, i_element) * ev_sqrt**power
      enddo
   enddo

   call dgemm('N','T',n_dim,n_dim,n_nonsingular,1._dp,work(:,1:n_nonsingular), &
               n_dim,work(:,1:n_nonsingular),n_dim,0._dp,matrix,n_dim)


   end subroutine power_genmat

  !> \brief Perform lapack diagonalization 
  !! @param[inout] matrix: working array
  !! @param[in] n_dim: dimension of the working array
  !! @param[in] threshold: threshold for the eigenvalue
  !! @param[in] n_singular: number of non singular eigenvalues
  !! @param[in] eigenvalues: array that contains the eigenvalues
  !! @param[in] work: auxiliary working array
   subroutine diagonalize_genmat(matrix,n_dim,threshold,n_nonsingular,eigenvalues,work) 

   integer                                 :: n_dim,n_nonsingular

   real(kind=dp)                            :: threshold
   real(kind=dp), dimension(n_dim)          :: eigenvalues
   real(kind=dp), dimension(n_dim,n_dim)    :: matrix, work

   ! Local variables

   integer                                 :: info, lwork
   integer, dimension(:), allocatable      :: ifail, iwork

   real(kind=dp)                            :: abs_tol, safe_minimum
   real(kind=dp), dimension(:), allocatable :: vwork

   ! Allocate work arrays
   allocate(vwork(8*n_dim),iwork(5*n_dim),ifail(n_dim))
   lwork = 8*n_dim

   ! Determine the safe minimum
   call get_machine_precision(safe_minimum)
   abs_tol = 2._dp*safe_minimum

   ! Perform diagonalization 
   call dsyevx ('V','V','U',n_dim,matrix,n_dim,threshold,1.d5,1,n_dim,abs_tol, &
                n_nonsingular,eigenvalues,work,n_dim,vwork,lwork,iwork,ifail,info)

   ! Print error messages 
   if (info.lt.0) then
      write(*,*) "Eigenvalue solver dsyevx:"
      write(*,*) "The ", -info,"th argument in dspevx had an illegal value."
   else if (info.gt.0) then
      write(*,*) "Eigenvalue solver dsyevx:", info,"eigenvectors failed to converge."
      stop
   end if

   deallocate(vwork,iwork,ifail) 
 
   end subroutine diagonalize_genmat


end module localized_basis_environments
