PROGRAM init_minimax

  USE kinds,       ONLY: dp
  USE init_types,  ONLY: input_type
  USE read_input,  ONLY: parse_command_line

  IMPLICIT NONE

  TYPE(input_type)                            :: input

  CALL parse_command_line(input)
 
  

END PROGRAM init_minimax
