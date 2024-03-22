program test_gx_localized_basis

   use kinds,                           only: dp        
   use localized_basis,                 only: gx_rirs_coefficients

   implicit none

   ! Local variables

   integer i_basbas, i_pair, i_point, unit_id
   integer n_basbas, n_basis_pairs, n_rk_points
   real(kind=dp)                              :: error
   real(kind=dp), dimension(:,:), allocatable :: ovlp2fn ! Testing
   real(kind=dp), dimension(:,:), allocatable :: ovlp3fn

   character(300) :: tmp_path, work_path

   !*******************************************************************!

   ! Initialization ! Hard coded

   n_basbas = 303
   n_basis_pairs = 1711
   n_rk_points = 638
   error = 0._dp

   ! Allocation of working arrays
   allocate(ovlp2fn(n_basis_pairs,n_rk_points))
   allocate(ovlp3fn(n_basis_pairs,n_basbas))

   ! Get temporary working path as an argument
   call get_command_argument(1, tmp_path)

   ! Get current working path
   call get_environment_variable('PWD',work_path)

   ! Read input data relative to the current working path 

   ! Read density computed on the real space grid
   open(newunit=unit_id,file=trim(adjustl(work_path))//'/test/localized_basis/density.dat')
   do i_pair=1,n_basis_pairs
      do i_point=1, n_rk_points
        read(unit_id,'(f28.18)') ovlp2fn(i_pair,i_point)
      end do
   end do
   close(unit=unit_id)

   ! Read three-center overlap coeeficients (RI-V)
   open(newunit=unit_id,file=trim(adjustl(work_path))//'/test/localized_basis/overlap.dat')
   do i_pair=1, n_basis_pairs
      do i_basbas=1, n_basbas
         read(unit_id,'(f28.18)') ovlp3fn(i_pair,i_basbas)
      end do
   end do
   close(unit_id)

   ! Call proccedure 
   call gx_rirs_coefficients(n_basbas,n_basis_pairs,n_rk_points, &
                             ovlp2fn,ovlp3fn,error)

   ! Deallocation of working arrays
   deallocate(ovlp2fn,ovlp3fn)

   ! Write the error into a temporary file to be read by the python driver 
   open(newunit=unit_id, file = trim(adjustl(tmp_path))//'/error.dat')
      write(unit_id,'(f30.15)') error
   close(unit_id)
   

end program test_gx_localized_basis
