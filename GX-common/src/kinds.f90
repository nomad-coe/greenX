!> Precision constants
module kinds
   implicit none
   public

   !> Single precision
   integer, parameter, public :: sp = selected_real_kind(6, 30)
   !> Double precision
   integer, parameter, public :: dp = selected_real_kind(14, 200)
   !> Medium length character
   integer, parameter, public :: medium_char = 100
   !> Long length character
   integer, parameter, public :: long_char = 200

   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: one = 1.0_dp

   character(len=1), parameter, public :: ch10 = char(10)

   public :: register_exc

! Private stuff
   integer,parameter :: err_len = 1024
   character(len=err_len), public, protected :: error_message__ = "No Error reported so far"

contains

subroutine register_exc(msg, filename, lineno)

!Arguments ------------------------------------
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: filename
 integer,optional,intent(in) :: lineno

! local variables
 integer :: my_lineno
 character(len=err_len) :: my_filename
! *********************************************************************

 my_filename = "Unknown File"; if (present(filename)) my_filename = trim(filename)
 my_lineno = 0; if (present(lineno)) my_lineno = lineno

 write(error_message__,'(7a,i0,7a)')ch10, &
   "--- ! ERROR", ch10, &
   "src_file: ", trim(my_filename), ch10, &
   "src_line: ",my_lineno, ch10, &
   "message: |",ch10,trim(msg),ch10, &
   "...",ch10

end subroutine register_exc

end module kinds
