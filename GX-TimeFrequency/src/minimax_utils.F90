module minimax_utils
#include "gx_common.h"
  use kinds, only: dp
  implicit none

  private

  public :: bsearch_erange

contains

  ! **************************************************************************************************
  !> \brief Modified bisection search to find first element in sorted array
  !>        that is strictly greater than a given value
  !>  lenght - lenght of sorted array
  !>  einter - sorted array of the energy intervals
  !>  eval - the energy value
  ! **************************************************************************************************
  function bsearch_erange(length, einter, eval) result(idx)
    integer, intent(in)                     :: length
    real(dp), dimension(length), intent(in) :: einter
    real(dp), intent(in)                    :: eval
    integer                                 :: idx

    ! Auxiliary variables
    integer                                 :: left, right, middle

    ! Begin work
    left = 1
    right = length
    idx = length + 1

    do while (left <= right)
       middle = (left + right) / 2
       if (einter(middle) <= eval) then
          left = middle + 1
       else
          right = middle - 1
          idx = middle
       end if
    end do

  end function bsearch_erange

end module minimax_utils
