module parse_minimax
    use kinds, only: dp, long_char
    use parse_utils, only: handle_io_stat
    implicit none
    
    private
    public :: read_minimax_eigenvalues, read_freq_grid

contains

  !> Read eigenvalues from file
  !>
  !> File must not have a header nor whitespace
  !>
  !> TODO  This routine needs to move
  !>       Extend to spin-polarized case 
  !>       Extend to n_kpoints
    subroutine read_minimax_eigenvalues(file_name, eigenvalues)

    !> Name of eigenvalue file
    character(len=*), intent(in) :: file_name
    !> Eigenvalues
    real(dp), intent(out) :: eigenvalues(:, :)

    integer :: i, stat

    if (size(eigenvalues, 2) /= 1) then
       write(*, *) 'read_minimax_eigenvalues not extended to spin-polarised case'
    endif

    open(unit=50, file=file_name,status='old',IOSTAT=stat)
    call handle_io_stat(stat, file_name)
    do i = 1, size(eigenvalues, 1)
       read(50, *) eigenvalues(i, 1)
    enddo
    close(50)

   end subroutine read_minimax_eigenvalues


  !> Read minimax frequency grid and weights
  subroutine read_freq_grid(npoints, file_name, erange, grid_freq, weights_freq,grid_time,weights_time)

    !> Number of minimax grid points
    integer, intent(in) :: npoints
    !> Minimax file name, prefixed by full path 
    character(len=*), intent(in) :: file_name
    !> Energy range (actually max energy diff / min energy diff)
    real(kind=dp), intent(in) :: erange
    
    !> Minimax grid for frequency and time
    real(kind=dp), dimension(:), intent(inout) :: grid_freq, grid_time
    !> Minimax weights for frequency and time
    real(kind=dp), dimension(:), intent(inout) :: weights_freq, weights_time

    character(len=long_char)                       :: line, ch 
    character(len=long_char), dimension(3)         :: tempStr
    integer                                        :: i, stat, istr
    integer                                        :: lastpos
    real(kind=dp)                                  :: lower_range, upper_range


    open(unit=62, file=trim(adjustl(file_name)), status='old', IOSTAT=stat)
    call handle_io_stat(stat, file_name)

    do
       read(62,'(A)',end=999) line                                ! read in whole line
       lastpos=-1                                                 ! start position when reading line from left to right
       istr = 0
       !*** find 'Erange' lines in file
       do i=1,len(line)
          ch=line(i:i)                                             ! read character at position i
          if(ch/='E'.and.i==1) EXIT                                ! if not start with E (i.e. Erange) skip
          if((ch==' '.or.i==len(line)).and.lastpos/= -1) then      ! if space: end of number or end of line?
             istr = istr + 1
             read(line(lastpos:i-1),*) tempStr(istr)                !lastpos is always start of string or number
             lastpos = -1
          endif
          if(ch/=' '.and.lastpos==-1) then
             lastpos=i
          end if
       enddo
       !*** read grid/weights if it is the right Erange
       if(istr > 0 ) then
          read(tempStr(2),*) lower_range
          read(tempStr(3),*) upper_range
          if(erange > lower_range.and.erange <= upper_range) then
             do i =  1, npoints
                if(i <= npoints) then
                   read(62,*) grid_freq(i), weights_freq(i), grid_time(i), weights_time(i)
                endif
             enddo
          endif
          EXIT
       endif

    enddo
999 continue
    close(62)

  end subroutine read_freq_grid

end module 
