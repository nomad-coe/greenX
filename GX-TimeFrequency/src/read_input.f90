MODULE read_input

 USE kinds,       ONLY: sp, dp,&
                        default_string_length
 USE init_types,  ONLY: input_type,&
                        minimax_type
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: parse_command_line, read_eigenvalues

  
CONTAINS

!****************************************************************
!> brief parsing of command line
!****************************************************************
SUBROUTINE parse_command_line(input)

  TYPE(input_type), INTENT(INOUT)                     :: input

  CHARACTER(LEN=default_string_length)                :: arg
  INTEGER(KIND=4)                                     :: narg
  INTEGER                                             :: i, stat

  input%have_grid_size = .FALSE.
  input%have_eigenvalues = .FALSE.
  input%have_nhomo = .FALSE.
  input%debug = .FALSE.
  input%output = .FALSE.

  narg=iargc()
  IF(narg < 1) THEN
      write(*,'(A76)') "Usage: ./init_minimax -s gridsize -e eigenvalue_file" //&
                       " -h nhomo -debug -output" 
      write(*,'(A38)') "Comment: -debug and -output are optional"
      STOP 
  ENDIF

  DO i=1, narg
     CALL getarg(i,arg)
     SELECT CASE(TRIM(arg))
     CASE("-s")
       input%have_grid_size = .TRUE.
       CALL getarg(i+1,input%npoints)
     CASE("-e")
       input%have_eigenvalues = .TRUE.
       CALL getarg(i+1,input%eigenvaluefile)
     CASE("-h")
       input%have_nhomo = .TRUE.
       CALL getarg(i+1,input%nhomo)
     CASE("-debug")
       input%debug = .TRUE.
     CASE("-output")
       input%output = .TRUE.
     END SELECT
  ENDDO

  IF(.NOT.(input%have_grid_size.AND.input%have_eigenvalues &
           .AND.input%have_nhomo)) THEN
     write(*,'(A62)') "Give the gridsize, file with DFT eigenvalue" //&
                      " and the HOMO state"
     STOP 
  ENDIF
 
END SUBROUTINE parse_command_line

!****************************************************************
!> brief reading eigenvalues from file, expects eigenvalues 
!>       in a.u. 
!****************************************************************
SUBROUTINE read_eigenvalues(input,minimax)

  TYPE(input_type), INTENT(IN)                     :: input
  TYPE(minimax_type), INTENT(INOUT)                :: minimax

  INTEGER                                          :: i, nlines, stat
   
  IF(input%have_nhomo) THEN
    READ(input%npoints,'(I6)') minimax%npoints
  ENDIF

  IF(input%have_nhomo) THEN
    READ(input%nhomo,'(I6)') minimax%nhomo
  ENDIF

  IF(input%have_eigenvalues)THEN
     nlines = 0
     OPEN(unit=50,file=input%eigenvaluefile,status='old',IOSTAT=stat)
     DO
       READ(50,*,iostat=stat)
       IF (stat/=0) EXIT
       nlines = nlines + 1
     END DO
     CLOSE(50)
     minimax%nstates =  nlines 

     ALLOCATE(minimax%eigenvalues(nlines))

     OPEN(unit=51,file=input%eigenvaluefile,status='old',IOSTAT=stat)
     DO i = 1, minimax%nstates
       READ(51,*) minimax%eigenvalues(i)
     ENDDO 
     CLOSE(51)
  ENDIF
 
END SUBROUTINE

END MODULE read_input 
