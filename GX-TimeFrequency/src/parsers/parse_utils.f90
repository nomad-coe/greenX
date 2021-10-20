!> Utility routines to assist in file parsing 
module parse_utils
    implicit none
    private 
    public :: get_n_lines, handle_io_stat

contains

 !> Get the number of lines in a file 
function get_n_lines(file_name) result(n_lines)
    !> File name
    character(len=*), intent(in) :: file_name
    !> Number of lines in file (including header and any white lines)
    integer :: n_lines

    !> Line counter
    integer :: line_counter
    !> IO status
    integer :: stat

    line_counter = 0
    open(unit=50, file=file_name, status='old', IOSTAT=stat)
    call handle_io_stat(stat, file_name)
    do
       read(50, *, iostat=stat)
       if (stat /= 0) exit
       line_counter = line_counter + 1
    enddo
    close(50)

    n_lines = line_counter

  end function get_n_lines

  !> Handle IO status 
  !>
  !> TODO(Alex) Add error handling
  subroutine handle_io_stat(stat, file_name)
    !> open status integer
    integer, intent(in) :: stat
    !> File name 
    character(len=*), intent(in) :: file_name
    if (stat /= 0) then
        write (*,*) "Open file:'"// trim(file_name) // "' failed with iostat = ", stat
        stop
    endif
  end subroutine handle_io_stat

end module parse_utils
