MODULE read_input

 USE kinds,       ONLY: sp, dp,&
                        default_string_length
 USE init_types,  ONLY: input_type,&
                        minimax_type
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: parse_command_line, read_eigenvalues,&
           read_freq_grid

  
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
  input%have_kpoints = .FALSE.
  input%have_gridpath = .FALSE.
  input%debug = .FALSE.
  input%output = .FALSE.

  narg=iargc()
  IF(narg < 1) THEN
      write(*,'(A98)') "Usage: ./init_minimax -s gridsize -e eigenvalue_file" //&
                       " -h nhomo -p path_to_grid_files -debug -output" 
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
     CASE("-p")
       input%have_gridpath = .TRUE.
       CALL getarg(i+1,input%filepath)
     CASE("-debug")
       input%debug = .TRUE.
     CASE("-output")
       input%output = .TRUE.
     END SELECT
  ENDDO

  IF(.NOT.(input%have_grid_size.AND.input%have_eigenvalues &
           .AND.input%have_nhomo.AND.input%have_gridpath)) THEN
     write(*,'(A90)') "Give the gridsize, file with DFT eigenvalue" //&
                      ", the HOMO state and the path to the grid files"
     STOP 
  ENDIF
 
END SUBROUTINE parse_command_line

!****************************************************************
!> brief reading eigenvalues from file, expects eigenvalues 
!>       in a.u. (assume spin-unpolarized case for now) 
!****************************************************************
SUBROUTINE read_eigenvalues(input,minimax)

  TYPE(input_type), INTENT(IN)                     :: input
  TYPE(minimax_type), INTENT(INOUT)                :: minimax

  INTEGER                                          :: i, nlines, stat
   
  IF(input%have_nhomo) THEN
    READ(input%npoints,'(I6)') minimax%npoints
  ENDIF

  IF(input%have_nhomo) THEN
    READ(input%nhomo,'(I6)') minimax%nhomo(1)
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

     ALLOCATE(minimax%eigenvalues(nlines,1))

     OPEN(unit=51,file=input%eigenvaluefile,status='old',IOSTAT=stat)
     DO i = 1, minimax%nstates
       READ(51,*) minimax%eigenvalues(i,1)
     ENDDO 
     CLOSE(51)
  ENDIF
 
END SUBROUTINE

!****************************************************************
!> brief reading frequency grids + weights
!****************************************************************
SUBROUTINE read_freq_grid(npoints,filepath,grid,weights,erange)

  INTEGER, INTENT(IN)                                        :: npoints
  CHARACTER(LEN=default_string_length), INTENT(IN)           :: filepath
  REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)                 :: grid
  REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)                 :: weights
  REAL(KIND=dp), INTENT(IN)                                  :: erange

  CHARACTER(LEN=default_string_length)                       :: filename
  CHARACTER(LEN=default_string_length)                       :: str_npoints
  INTEGER                                                    :: i, stat

  WRITE(str_npoints,*) npoints

  filename = TRIM(filepath) // TRIM(ADJUSTL(str_npoints)) // "_freq_points.dat"
  OPEN(unit=62,file=filename,status='old',IOSTAT=stat)
  ! now start reading the file at the right position

  !DO i =  1, 2*npoints
  !  IF(i.le.npoints) THEN
  !    READ(62,*) grid(i)
  !  ELSE
  !    READ(62,*) weights(i)
  !  ENDDO
  !ENDDO
  CLOSE(62)

END SUBROUTINE read_freq_grid

END MODULE read_input 
