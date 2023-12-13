! ***************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

module gx_ac
   use kinds, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_double_complex, c_ptr
   use pade_approximant, only: evaluate_thiele_pade, thiele_pade, c_zero, c_one
   implicit none

   public :: thiele_pade_api, &
      params_mp, &
      create_thiele_pade_mp, &
      evaluate_thiele_pade_mp, &
      free_pade_model

   !> brief store arbitrary precision parameters
   type :: params_mp
      logical     :: initialized = .False.
      type(c_ptr) :: params_ptr
   end type params_mp

   interface

      !> auxiliary function to compute Thiele-Pade parameters using arbitrary precision numbers
      !! @param[in] n_par - order of the interpolant
      !! @param[in] x_ref - array of the reference points
      !! @param[in] y_ref - array of the reference function values
      !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
      !! @return - pointer to abstract type to store all parameters
      function thiele_pade_mp_aux(n_par, x_ref, y_ref, do_greedy) bind(C, name="thiele_pade_mp")
         import :: c_double_complex, c_int, c_ptr
         integer(c_int), value                    :: n_par
         complex(c_double_complex), dimension(*)  :: x_ref
         complex(c_double_complex), dimension(*)  :: y_ref
         integer(c_int), value                    :: do_greedy
         type(c_ptr)                              :: thiele_pade_mp_aux
      end function thiele_pade_mp_aux

      !> auxiliary function to evaluate the Thiele-Pade parameters using arbitrary precision numbers
      !! @param[in] x - point where the function needs to be evaluated
      !! @param[in] params_ptr - pointer to abstract type to store all parameters
      !! @return - interpolated function value
      function evaluate_thiele_pade_mp_aux(x, params_ptr) bind(C, name="evaluate_thiele_pade_mp")
         import :: c_ptr, c_double_complex
         complex(c_double_complex), value    :: x
         type(c_ptr), value                  :: params_ptr
         complex(c_double_complex)           :: evaluate_thiele_pade_mp_aux
      end function evaluate_thiele_pade_mp_aux

      !> Frees the C++ pointers used in the Fortran type
      !! @param[in] params_ptr - the C++ pointer
      subroutine free_pade_model(params_ptr) bind(C, name="free_pade_model")
         import :: c_ptr
         type(c_ptr), value :: params_ptr
      end subroutine free_pade_model

   end interface

contains

   !> API function to compute Thiele-Pade approximations of a meromorphic function
   !! @param[in] n_par - order of the interpolant
   !! @param[in] x_ref - array of the reference points
   !! @param[in] y_ref - array of the reference function values
   !! @param[in] x_query - array of points where the function needs to be evaluated
   !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
   !! @param[out] y_query - array of the interpolated values at x_query
   subroutine thiele_pade_api(n_par, x_ref, y_ref, x_query, y_query, do_greedy)
      integer, intent(in)                         :: n_par
      complex(kind=dp), dimension(:), intent(in)  :: x_ref, y_ref, x_query
      complex(kind=dp), dimension(:), intent(out) :: y_query
      logical, optional, intent(in)               :: do_greedy

      ! Internal variables
      integer                                     :: i, num_query
      complex(kind=dp), dimension(size(x_ref))    :: x_ref_local
      complex(kind=dp), dimension(n_par)          :: a_par

      ! Compute the coefficients a_par
      x_ref_local(:) = x_ref
      call thiele_pade(n_par, x_ref_local, y_ref, a_par, do_greedy)

      ! Compute the number of query points
      num_query = size(x_query)

      ! Evaluate the Thiele-Pade approximation at the query points
      do i = 1, num_query
         call evaluate_thiele_pade(n_par, x_ref_local, x_query(i), a_par, y_query(i))
      end do

   end subroutine thiele_pade_api

   !> API function to compute Thiele-Pade parameters using arbitrary precision numbers
   !! @param[in] n_par - order of the interpolant
   !! @param[in] x_ref - array of the reference points
   !! @param[in] y_ref - array of the reference function values
   !! @param[out] params - abstract type to store all parameters in arb. prec. representation
   !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
   subroutine create_thiele_pade_mp(n_par, x_ref, y_ref, params, do_greedy)
      integer, intent(in)                             :: n_par
      complex(kind=dp), dimension(:), intent(in)      :: x_ref, y_ref
      type(params_mp), intent(out)                    :: params
      logical, optional, intent(in)                   :: do_greedy

      ! Internal variables
      integer                                         :: local_do_greedy = 1

      ! use integer bools for interoperability with C
      if (present(do_greedy)) then
         if (do_greedy) then
            local_do_greedy = 1
         else
            local_do_greedy = 0
         end if
      end if

      ! compute coefficients
      params%initialized = .True.
      params%params_ptr = thiele_pade_mp_aux(n_par, x_ref, y_ref, local_do_greedy)

   end subroutine create_thiele_pade_mp

   !> API function to evaluate the Thiele-Pade approximation using arbitrary precision numbers
   !! @param[in] x - point where the function is evaluated
   !! @param[out] y - interpolated function value at x
   !! @param[in] params - abstract type to store all parameters in arb. prec. representation
   subroutine evaluate_thiele_pade_mp(x, y, params)
      complex(kind=dp), intent(in)     :: x
      complex(kind=dp), intent(out)    :: y
      type(params_mp), intent(in)      :: params
      
      y = cmplx(0.0d0, 0.0d0, kind=8)
      if (.not. params%initialized)  return

      y = evaluate_thiele_pade_mp_aux(x, params%params_ptr)

   end subroutine evaluate_thiele_pade_mp

end module gx_ac
