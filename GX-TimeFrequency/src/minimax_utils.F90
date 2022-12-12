! ***************************************************************************************************
!  Copyright (C) 2020-2022 Green-X library                                                          
!  This file is distributed under the terms of the APACHE2 License.                                 
!                                                                                                   
! ***************************************************************************************************
!> \brief This module contains auxiliary procedures and data structures for the main minimax routines
! ***************************************************************************************************
module minimax_utils
#include "gx_common.h"
  use kinds, only: dp
  implicit none

  private

  type :: er_aw_aux
     ! Sorted array of the energy ranges
     real(kind=dp), dimension(:), allocatable :: energy_range
     ! Matrices with coefficients and weights per energy region
     real(kind=dp), dimension(:, :), allocatable :: aw_erange_matrix
  end type er_aw_aux

  public :: bsearch_erange, er_aw_aux

contains

  !> \brief Modified bisection search to find first element in sorted array
  !>        that is strictly greater than a given value
  !! @param[in] lenght - lenght of sorted array
  !! @param[in] einter - sorted array of the energy intervals
  !! @param[in] eval - the energy value
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
