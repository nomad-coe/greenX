MODULE init_types
 
 USE kinds,    ONLY: dp, default_string_length
 
 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: input_type, minimax_type

!***************************************************************************
 TYPE input_type
   LOGICAL                                  :: have_grid_size
   LOGICAL                                  :: have_eigenvalues
   LOGICAL                                  :: have_nhomo
   LOGICAL                                  :: debug
   LOGICAL                                  :: output
   CHARACTER(LEN=default_string_length)     :: eigenvaluefile
   CHARACTER(LEN=default_string_length)     :: nhomo, npoints
 END TYPE input_type
!***************************************************************************
 TYPE minimax_type
   INTEGER                                             :: nhomo
   INTEGER                                             :: npoints
   INTEGER                                             :: nstates
   REAL(KIND=dp)                                       :: erange
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: eigenvalues
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: grid
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE            :: weights 
 END TYPE minimax_type

!***************************************************************************

END MODULE init_types
