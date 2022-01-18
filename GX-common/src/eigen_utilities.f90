module eigen_utilities
  use kinds, only: sp, dp
  implicit none
  private
  public :: get_erange_non_periodic

contains

  !> Calculate the energy range for non-periodic systems
  !>
  !> emin is defined as min(Eg_spinup, Eg_spindown).
  !> emax is defined as max([E_max - E_min]_spinup, [E_max - E_min]_spindown)
  !>
  !> TODO(Alex) Why is this called erange when a ratio is returned?
  !> Rename emin as erange_min, and emax as erange_max
  function get_erange_non_periodic(eigenvalues, nhomo) result(erange)

    !> Eigenvalues
    real(dp), dimension(:, :), intent(in) :: eigenvalues
    !> HOMO index
    integer, dimension(:), intent(in) :: nhomo
    !> Energy range
    real(dp) :: erange
    !> Arbitrarily large value
    real(dp), parameter :: large_energy = 1.E50_dp

    integer :: ispin, nspin, ihomo, ilumo
    real(dp) :: emin, emax, emin_old, emax_old

    nspin = size(eigenvalues, 2)
    emax_old = -large_energy
    emin_old = large_energy

    do ispin = 1, nspin
      ihomo = nhomo(ispin)
      ilumo = nhomo(ispin) + 1
      emin = eigenvalues(ilumo, ispin) - eigenvalues(ihomo, ispin)
      emax = maxval(eigenvalues(:, ispin)) - minval(eigenvalues(:, ispin))
      emin = min(emin, emin_old)
      emax = max(emax, emax_old)
      emin_old = emin
      emax_old = emax
    enddo
    erange = emax/emin

  end function get_erange_non_periodic

end module eigen_utilities
