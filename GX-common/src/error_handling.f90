module error_handling
   implicit none
   private

   !> Arbitrary error length
   integer, parameter :: err_len = 1024
   !> Error message
   character(len=err_len), public, protected :: error_message__ = "No Error reported so far"
   !> TODO Move me
   character(len=1), parameter, public :: ch10 = char(10)

contains

   !> Write an error to file.
   !>
   !> Note, this should replaced with somethign machine-readable.
   subroutine register_exc(msg, filename, line_number)

      !>  Error message
      character(len=*), intent(in) :: msg
      !> File name
      character(len=*), optional, intent(in) :: filename
      !> Line number
      integer, optional, intent(in) :: line_number

      integer :: line_num
      character(len=err_len) :: my_filename

      my_filename = "Unknown File"
      if (present(filename)) then
         my_filename = trim(filename)
      end if

      line_num = 0
      if (present(line_number)) then
         line_num = line_number
      end if

      write (error_message__, '(7a,i0,7a)') ch10, &
         "--- ! ERROR", ch10, &
         "src_file: ", trim(my_filename), ch10, &
         "src_line: ", line_num, ch10, &
         "message: |", ch10, trim(msg), ch10, &
         "...", ch10

   end subroutine register_exc

end module error_handling
