MODULE read_input

 USE kinds,       ONLY: sp, dp,&
                        default_string_length
 USE init_types,  ONLY: input_type
                     ! matrix_structure_type
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: parse_command_line!, read_or_generate_matrix

  
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

   IF(.NOT.(input%have_grid_size.OR.input%have_eigenvalues .OR.input%have_nhomo)) THEN
      write(*,'(A50)') "Give the gridsize, File with DFT eigenvalue and the HOMO state"
      STOP 
   ENDIF
  ENDDO
 
END SUBROUTINE parse_command_line

!SUBROUTINE read_r_generate_matrix(input,structure,matrix_A,matrix_B)
!
! TYPE(input_type), INTENT(IN)                     :: input
! TYPE(matrix_structure_type), INTENT(INOUT)       :: structure
! REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE,&
!   INTENT(INOUT)                                  :: matrix_A, matrix_B   
!
! INTEGER                                          :: rowA, colA,&
!                                                     rowB, colB
! INTEGER                                          :: irow, icol, stat
! INTEGER                                          :: seed
!
! IF(input%have_matrix_file)THEN
!    OPEN(unit=50,file=input%matrix_file,status='old',IOSTAT=stat)
!
!
!    READ(50,*) rowA, colA
!    
!    ALLOCATE(matrix_A(rowA,colA))
!
!    DO icol=1,colA
!     DO irow=1,rowA
!        READ(50,*) matrix_A(irow,icol)
!     END DO
!    END DO
!
!    READ(50,*) rowB, colB
!
!    ALLOCATE(matrix_B(rowB,colB))
!
!    IF(colA /= rowB) THEN
!      WRITE(*,*)'column A not equal row B'
!      STOP
!    ENDIF
!
!    DO icol=1,colB
!     DO irow=1,rowB
!        READ(50,*) matrix_B(irow,icol)
!     END DO
!    END DO
!
!    CLOSE(50)
! ELSE 
!    READ(input%mdimension,'(I6)')  rowA
!    READ(input%mdimension,'(I6)')  colA
!    READ(input%mdimension,'(I6)')  rowB
!    READ(input%mdimension,'(I6)')  colB
!    ALLOCATE(matrix_A(rowA,colA))
!    ALLOCATE(matrix_B(rowB,colB))
!    seed=1299
!    DO icol=1,colA
!     DO irow=1,rowA
!        matrix_A(irow,icol) = RAND(seed)
!        seed=seed+1000*irow*icol
!     END DO
!    END DO
!
!    DO icol=1,colB
!     DO irow=1,rowB
!        matrix_B(irow,icol) = RAND(seed)
!        seed=seed+1000*irow*icol
!     END DO
!    END DO
! ENDIF
!
!  structure%rowA=rowA 
!  structure%colA=colA 
!  structure%rowB=rowB 
!  structure%colB=colB
!
!  IF(input%have_block) THEN
!    READ(input%mblock,'(I6)')  structure%nblock
!  ELSE
!    structure%nblock = 16
!  ENDIF 

!END SUBROUTINE

END MODULE read_input
