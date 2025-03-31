! **************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief This module contains the subroutines for the localized basis set component of the library
! ***************************************************************************************************
module localized_basis_environments

   use kinds,                 only: dp
   use error_handling,        only: register_exc   
   use localized_basis_types, only: separable_ri_types, &
                                    polarizability_types
   use lapack_interfaces,     only: dgemm, dsyevx

   implicit none

   !*******************************************************************!  

   contains

  !> \brief Initialize the ri_rs types.
  !! @param[inout] ri_rs: type for the separable ri
   subroutine initialization(ri_rs)

   type(separable_ri_types) :: ri_rs

   ri_rs%basis%n_basbas      =303
   ri_rs%basis%n_basis       = 58 
   ri_rs%basis%n_loc_basbas  = 303
   ri_rs%basis%n_basis_pairs = 1711
   ri_rs%n_points            = 638
   ri_rs%error               = 0.d0

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

  ! **********************************************************************
  !> \brief Initialize kohn sham array
  !  param[inout] pi_pq:  polarizability environment
  ! **********************************************************************   
   subroutine initialize_kohn_sham (pi_pq)

     type(polarizability_types)               :: pi_pq

     if (.not.allocated(pi_pq%ks%eigenvalues)) then
        allocate(pi_pq%ks%eigenvalues(pi_pq%ks%n_states, 1))
     end if

     if (.not.allocated(pi_pq%ks%eigenvectors)) then
        allocate(pi_pq%ks%eigenvectors(pi_pq%ks%n_basis, pi_pq%ks%n_states, 1))
     end if

     if (.not.allocated(pi_pq%ks%wave)) then
        allocate(pi_pq%ks%wave(pi_pq%ks%n_basis, pi_pq%ri_rs%n_points, 1))
     end if

   end subroutine initialize_kohn_sham        

  ! *********************************************************************
  !> brief allocate the minimax time-frequency grids 
  !  param[inout] pi_pq:  polarizability environment
  ! *********************************************************************
   subroutine initialize_minimax_grids(pi_pq)

     type(polarizability_types)               :: pi_pq

     !Local variables
     integer                                  :: n_points

     ! Hard coded (change later)
     n_points = 6

     if (.not.allocated(pi_pq%minimax%cos_tf)) then
        allocate(pi_pq%minimax%cos_tf(n_points, n_points))
     end if

     if (.not.allocated(pi_pq%minimax%omega)) then
        allocate(pi_pq%minimax%omega(n_points))
     end if

     if (.not.allocated(pi_pq%minimax%tau)) then
        allocate(pi_pq%minimax%tau(n_points))
     end if

     if (.not.allocated(pi_pq%minimax%weights)) then
        allocate(pi_pq%minimax%weights(n_points))
     end if

     pi_pq%minimax%n_points    = n_points
     pi_pq%minimax%cos_tf(:,:) = 0.d0
     pi_pq%minimax%omega(:)    = 0.d0
     pi_pq%minimax%tau(:)      = 0.d0
     pi_pq%minimax%weights(:)  = 0.d0
     

   end subroutine initialize_minimax_grids

  ! **********************************************************************
  !> \brief Initalize polarizability array   
  !  param[inout] pi_pq:  polarizability environment
  ! **********************************************************************   
   subroutine initialize_polarizability(pi_pq)

     type(polarizability_types)               :: pi_pq

     ! Initialize polarizability arrays
     if (.not.allocated(pi_pq%tau)) then
        allocate(pi_pq%tau(pi_pq%ri_rs%basis%n_basbas,pi_pq%ri_rs%basis%n_loc_basbas))
     end if

     if (.not.allocated(pi_pq%omega)) then
        allocate(pi_pq%omega(pi_pq%ri_rs%basis%n_basbas,pi_pq%ri_rs%basis%n_loc_basbas, &
                             pi_pq%minimax%n_points))
     end if

     if (.not.allocated(pi_pq%chi%matrix)) then
        allocate(pi_pq%chi%matrix(pi_pq%ri_rs%n_points,pi_pq%ri_rs%n_points))
     end if

     pi_pq%chi%matrix(:,:) = 0.d0
     pi_pq%omega(:,:,:)    = 0.d0
     pi_pq%tau(:,:)        = 0.d0  

   end subroutine initialize_polarizability        
     
  !> \brief Deallocate the ri_rs types.
  !! @param[inout] ri_rs: type for the separable ri   
   subroutine deallocations(ri_rs, keep_coeff)

   type(separable_ri_types) :: ri_rs

   logical, optional        :: keep_coeff

   ! Local variables 
   logical                  ::  my_keep_coeff

   if (present(keep_coeff)) then
      my_keep_coeff = keep_coeff
   else
      my_keep_coeff = .false.
   end if

   if (allocated(ri_rs%ovlp2fn)) deallocate(ri_rs%ovlp2fn)

   if (allocated(ri_rs%ovlp3fn)) deallocate(ri_rs%ovlp3fn)

   if (.not.my_keep_coeff) then 
      if (allocated(ri_rs%z_coeff)) deallocate(ri_rs%z_coeff)
   end if

   end subroutine deallocations

  ! **********************************************************************
  !> \brief Deallocate minimax array   
  !  param[inout] pi_pq:  polarizability environment
  ! **********************************************************************
   subroutine deallocate_minimax_grids(pi_pq)

     type(polarizability_types)               :: pi_pq

     if (allocated(pi_pq%minimax%cos_tf))  deallocate(pi_pq%minimax%cos_tf)

     if (allocated(pi_pq%minimax%omega))    deallocate(pi_pq%minimax%omega)

     if (allocated(pi_pq%minimax%tau))     deallocate(pi_pq%minimax%tau)
     
     if (allocated(pi_pq%minimax%weights)) deallocate(pi_pq%minimax%weights)

   end subroutine deallocate_minimax_grids

    ! **********************************************************************
  !> \brief Deallocate kohn sham arrays
  !  param[inout] pi_pq:  polarizability environment
  ! **********************************************************************   
   subroutine deallocate_kohn_sham (pi_pq)

     type(polarizability_types)               :: pi_pq

     if (allocated(pi_pq%ks%eigenvalues))  deallocate(pi_pq%ks%eigenvalues)

     if (allocated(pi_pq%ks%eigenvectors)) deallocate(pi_pq%ks%eigenvectors)

     if (allocated(pi_pq%ks%wave))          deallocate(pi_pq%ks%wave)

   end subroutine deallocate_kohn_sham

  ! **********************************************************************
  !> \brief Initalize polarizability array   
  !  param[inout] pi_pq:  polarizability environment
  ! **********************************************************************   
   subroutine deallocate_polarizability(pi_pq, keep_pi)

     type(polarizability_types) :: pi_pq

     logical, optional          :: keep_pi

     ! Local variables 
     logical                    ::  my_keep_pi

     if (present(keep_pi)) then
        my_keep_pi = keep_pi
     else
        my_keep_pi = .false.
     end if

     call deallocate_minimax_grids(pi_pq)

     call deallocate_kohn_sham (pi_pq)

     if (allocated(pi_pq%tau))        deallocate(pi_pq%tau)

     if(.not.my_keep_pi) then
       if (allocated(pi_pq%omega))      deallocate(pi_pq%omega)
     end if

     if (allocated(pi_pq%chi%matrix)) deallocate(pi_pq%chi%matrix)

   end subroutine deallocate_polarizability   

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
