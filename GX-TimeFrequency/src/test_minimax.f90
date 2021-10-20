!> Test Program, demonstrating how to call minimax grid routine
!> 
!> Example run command: ./exe -s 32 -e ../../GX-TimeFrequency/src/ref/H2O_eigenvalues.txt -h 5 -output 32_freq_points.dat
!> TODO
!> * Implement grid scaling
!> * Write a library module to expose the object: 'use GXTimeFrequency, only: minimax' 
program test_minimax
  use parse_minimax,   only: read_minimax_eigenvalues !! "Some" eigenvalue file parser
  use parse_utils,     only: get_n_lines !! Utility to get number of lines in a file
  use read_cmd_line,   only: input_type  !! Makes running the test for multiple files easier
  use minimax,         only: minimax_type !! Object under test
  implicit none

  !> Command line options
  type(input_type) :: input
  !> Minimax grid
  type(minimax_type) :: mm
  !> Number of grid points, as a string 
  character(len=100) :: str_npoints
  
  call input%parse()
  call mm%init(input%npoints, get_n_lines(input%eigenvaluefile), input%nhomo)
  call read_minimax_eigenvalues(input%eigenvaluefile, mm%eigenvalues)
  call mm%compute_e_range(input%have_kpoints)
  call mm%read_frequency_weights()
  if (input%output) then
    write(str_npoints, *) mm%n_points
    call mm%write(trim(adjustl(input%output_file)))
  endif
  call mm%deallocate()

end program test_minimax
