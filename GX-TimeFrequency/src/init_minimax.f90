PROGRAM init_minimax

  USE kinds,       ONLY: dp
  USE init_types,  ONLY: input_type,&
                         minimax_type
  USE read_input,  ONLY: parse_command_line,&
                         read_eigenvalues

  IMPLICIT NONE

  TYPE(input_type)                            :: input
  TYPE(minimax_type)                          :: minimax

  CALL parse_command_line(input)
  CALL read_eigenvalues(input,minimax)
 
  

END PROGRAM init_minimax
