module localized_basis_environments

   use kinds,             only: dp
   use localized_basis_types
   use lapack_interfaces, only: dgemm, dsyevx


   implicit none

   !*******************************************************************!  

   contains

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

   subroutine deallocations(ri_rs)

   type(separable_ri_types) :: ri_rs

   if (allocated(ri_rs%ovlp2fn)) deallocate(ri_rs%ovlp2fn)

   if (allocated(ri_rs%ovlp3fn)) deallocate(ri_rs%ovlp3fn)

   if (allocated(ri_rs%z_coeff)) deallocate(ri_rs%z_coeff)   

   end subroutine deallocations

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

   subroutine get_machine_precision(safe_minimum)

   real(kind=dp) :: safe_minimum
   real(kind=dp), external ::  dlamch

   safe_minimum = dlamch ('S') 

   end subroutine get_machine_precision

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
