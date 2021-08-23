PROGRAM init_minimax

  USE kinds,       ONLY: dp
  USE init_types,  ONLY: input_type,&
                         minimax_type,&
                         allocate_minimax_structure,&
                         deallocate_minimax_structure
  USE read_input,  ONLY: parse_command_line,&
                         read_eigenvalues,&
                         read_freq_grid
  USE init_grids,  ONLY: get_erange_non_periodic

  IMPLICIT NONE

  TYPE(input_type)                            :: input
  TYPE(minimax_type)                          :: minimax
  REAL(KIND=dp)                               :: erange

  CALL parse_command_line(input)
  CALL read_eigenvalues(input,minimax)
  CALL allocate_minimax_structure(minimax)
  IF(input%have_kpoints) THEN
    ! needs to be still added
  ELSE
    erange = get_erange_non_periodic(minimax%eigenvalues,minimax%nhomo)
  ENDIF
  CALL read_freq_grid(minimax%npoints,input%filepath,minimax%grid,&
                      minimax%weights,erange)
  
  !CALL scale_grids...
  CALL deallocate_minimax_structure(minimax)

END PROGRAM init_minimax
