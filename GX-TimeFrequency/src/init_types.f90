MODULE init_types
 
 USE kinds,    ONLY: dp, default_string_length
 
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: input_type, minimax_type, allocate_minimax_structure,&
           deallocate_minimax_structure

!***************************************************************************
 TYPE input_type
   LOGICAL                                  :: have_grid_size
   LOGICAL                                  :: have_eigenvalues
   LOGICAL                                  :: have_nhomo
   LOGICAL                                  :: have_kpoints
   LOGICAL                                  :: have_gridpath
   LOGICAL                                  :: debug
   LOGICAL                                  :: output
   CHARACTER(LEN=default_string_length)     :: filepath
   CHARACTER(LEN=default_string_length)     :: eigenvaluefile
   CHARACTER(LEN=default_string_length)     :: nhomo, npoints
 END TYPE input_type
!***************************************************************************
 TYPE minimax_type
   INTEGER                                             :: npoints
   INTEGER                                             :: nstates
   INTEGER, DIMENSION(2)                               :: nhomo
   REAL(KIND=dp)                                       :: erange
   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE          :: eigenvalues
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: grid
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: weights 
 END TYPE minimax_type
!***************************************************************************

CONTAINS

!***************************************************************************
!> brief allocate minimax array
!***************************************************************************
SUBROUTINE allocate_minimax_structure(minimax)

  TYPE(minimax_type)                          :: minimax

  ALLOCATE(minimax%grid(minimax%npoints)) 
  ALLOCATE(minimax%weights(minimax%npoints))

  minimax%grid(:) = 0.0_dp 
  minimax%weights(:) = 0.0_dp 

END SUBROUTINE allocate_minimax_structure

!***************************************************************************
!> brief deallocate minimax array
!***************************************************************************
SUBROUTINE deallocate_minimax_structure(minimax)

  TYPE(minimax_type)                          :: minimax

  IF(ALLOCATED(minimax%eigenvalues)) DEALLOCATE(minimax%eigenvalues)
  IF(ALLOCATED(minimax%grid)) DEALLOCATE(minimax%grid)
  IF(ALLOCATED(minimax%weights)) DEALLOCATE(minimax%weights)

END SUBROUTINE deallocate_minimax_structure

END MODULE init_types
