MODULE init_grids

 USE kinds,       ONLY: sp, dp,&
                        default_string_length

 IMPLICIT NONE
 
 PRIVATE

 PUBLIC :: get_erange_non_periodic

  
CONTAINS

!*********************************************************************************
!> calculate the energy range for non-periodic systems
!*********************************************************************************
FUNCTION get_erange_non_periodic(eigenvalues,nhomo) RESULT(erange)

  REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)           :: eigenvalues
  INTEGER, DIMENSION(:), INTENT(IN)                   :: nhomo
  REAL(KIND=dp)                                       :: erange

  INTEGER                                             :: ispin, nspin
  REAL(KIND=dp)                                       :: emin, emax,&
                                                         emin_old, emax_old


  nspin = SIZE(eigenvalues,2)
  emax_old = -1.E50_dp
  emin_old = 1.E50_dp

  DO ispin = 1, nspin
    emin = eigenvalues(nhomo(ispin) + 1, ispin) - eigenvalues(nhomo(ispin),ispin)
    emax = MAXVAL(eigenvalues(:,ispin)) -  MINVAL(eigenvalues(:,ispin))
    emin = MIN(emin,emin_old)
    emax = MAX(emax,emax_old)
    emin_old = emin
    emax_old = emax
  ENDDO
  erange = emax/emin 


END FUNCTION get_erange_non_periodic
                
END MODULE init_grids
