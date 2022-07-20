!> Precision constants
module kinds
   implicit none
   private

   !> Single precison
   integer, parameter, public :: sp = selected_real_kind(6, 30)
   !> Double precision
   integer, parameter, public :: dp = selected_real_kind(14, 200)
   !> Medium length character
   integer, parameter, public :: medium_char = 100
   !> Long length character
   integer, parameter, public :: long_char = 200

   character(len=1), parameter :: ch10 = char(10)

   public :: register_exc

! Private stuff
   character(len=4*1024) :: error_message__ = "No Error reported so far"

contains

subroutine register_exc(msg, filename, lineno)

!Arguments ------------------------------------
  character(len=*), intent(in) :: msg, filename
  integer,optional,intent(in) :: lineno
! *********************************************************************

  write(error_message__,'(7a,i0,7a)')ch10, &
    "--- ! ERROR", ch10, &
    "src_file: ", trim(filename), ch10, &
    "src_line: ",lineno, ch10, &
    "message: |",ch10,trim(msg),ch10, &
    "...",ch10

end subroutine register_exc

end module kinds
