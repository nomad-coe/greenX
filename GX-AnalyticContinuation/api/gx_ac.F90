! ***************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

module gx_ac
   use kinds, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_double_complex, c_ptr
   use pade_approximant, only: evaluate_thiele_pade, thiele_pade, c_zero, c_one
   implicit none

   private 
   public :: thiele_pade_api, &
             params, &
             create_thiele_pade, &
             evaluate_thiele_pade_at, &
             free_params, &
             arbitrary_precision_available

   !> brief store the parameters of the tiehle pade model 
   !> (potentially in abitrary precision floats using GMP)
   type :: params 
       logical    :: initialized = .false.
       integer    :: n_par
       integer    :: precision 
       logical    :: multiprecision_used_internally

       ! for multiple precision arithmetic
       type(c_ptr) :: params_ptr

       ! for double precision arithmetic
       complex(kind=dp), dimension(:), allocatable :: x_ref
       complex(kind=dp), dimension(:), allocatable :: a_par
   end type params 

#ifdef GMPXX_FOUND
   logical, parameter :: arbitrary_precision_available = .true.
#else
   logical, parameter :: arbitrary_precision_available = .false.
#endif

#ifdef GMPXX_FOUND
   interface

      !> auxiliary function to compute Thiele-Pade parameters using arbitrary precision numbers
      !! @param[in] n_par - order of the interpolant
      !! @param[in] x_ref - array of the reference points
      !! @param[in] y_ref - array of the reference function values
      !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
      !! @return - pointer to abstract type to store all parameters
      function thiele_pade_mp_aux(n_par, x_ref, y_ref, do_greedy, precision) bind(C, name="thiele_pade_mp")
         import :: c_double_complex, c_int, c_ptr
         integer(c_int), value                    :: n_par
         complex(c_double_complex), dimension(*)  :: x_ref
         complex(c_double_complex), dimension(*)  :: y_ref
         integer(c_int), value                    :: do_greedy
         integer(c_int), value                    :: precision
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
#endif

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



   !> API function to compute Thiele-Pade parameters 
   !! (potentially using arbitrary precision arithmetic)
   !!
   !! @param[in] n_par - order of the interpolant
   !! @param[in] x     - array of the reference points
   !! @param[in] y     - array of the reference function values
   !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
   !! @param[in] precision - precision in bits (!! not bytes !!)
   !! @return    params - abstract type to store all parameters in arb. prec. representation
   type(params) function create_thiele_pade(n_par, x, y, do_greedy, precision) result(par)
      integer, intent(in)                        :: n_par
      complex(kind=dp), dimension(:), intent(in) :: x, y
      logical, optional, intent(in)              :: do_greedy
      integer, optional, intent(in)              :: precision 

      ! Internal variables
      integer                                         :: local_do_greedy = 1

      ! initialize type
      par%initialized = .true.
      par%n_par = n_par

      ! precision of arithmetic internally
      if (present(precision)) then 
        if (precision .eq. 64) then 
            ! double precision case
            par%precision = precision
            par%multiprecision_used_internally = .false.
        else 
            ! arbitrary precision case
            par%precision = precision
            par%multiprecision_used_internally = .true.
#ifdef GMPXX_FOUND
#else 
            print *, "*** ERROR: multiple precision float arithmetic requested but not linked against GMP library"
            return 
#endif
            
        end if 
      else 
#ifdef GMPXX_FOUND
        ! default is quadrupel precision if GMP is linked
        par%precision = 128
        par%multiprecision_used_internally = .true.
#else 
        ! default is double precision if GMP is not linked
        par%precision = 64
        par%multiprecision_used_internally = .false.
#endif
      end if 

      if (.not. par%multiprecision_used_internally) then 

        allocate(par%a_par(n_par))
        allocate(par%x_ref(size(x)))

        ! compute the coefficients 
        par%x_ref(:) = x
        call thiele_pade(n_par, par%x_ref, y, par%a_par, do_greedy)

      elseif (par%multiprecision_used_internally) then 

        ! use integer bools for interoperability with C
        if (present(do_greedy)) then
            if (do_greedy) then
                local_do_greedy = 1
            else
                local_do_greedy = 0
            end if
        end if
#ifdef GMPXX_FOUND
        par%params_ptr = thiele_pade_mp_aux(n_par, x, y, local_do_greedy, par%precision)
#endif

      end if 

   end function create_thiele_pade 



   !> API function to evaluate the Thiele-Pade approximation 
   !! (potentially using arbitrary precision numbers)
   !!
   !! @param[in] x - point where the function is evaluated
   !! @param[in] params - abstract type to store all parameters 
   !! @return    y - interpolated function value at x
   function evaluate_thiele_pade_at(par, x) result(y)
      type(params), intent(in) :: par
      complex(kind=dp), dimension(:), intent(in) :: x
      complex(kind=dp), dimension(size(x)) :: y

      ! internal variables
      integer :: num_query, i

      ! initialized?
      if (par%initialized .eqv. .false.) then 
          print *, "WARNING: pade parameters not initialized"
          y(:) = 0.0d0 
          return 
      end if 

      ! Compute the number of query points
      num_query = size(x)

      ! Evaluate the Thiele-Pade approximation at the query points
      do i = 1, num_query

        if (.not. par%multiprecision_used_internally) then 
            call evaluate_thiele_pade(par%n_par, par%x_ref, x(i), par%a_par, y(i))
        elseif (par%multiprecision_used_internally) then 
#ifdef GMPXX_FOUND
            y(i) = evaluate_thiele_pade_mp_aux(x(i), par%params_ptr)
#endif
        end if 

      end do

   end function evaluate_thiele_pade_at



   !> deallocate the pade model 
   !!
   !! @param[inout] par - pade model struct 
   subroutine free_params(par)
       type(params), intent(inout) :: par 

       if (.not. par%initialized) return 

       if (allocated(par%a_par)) deallocate(par%a_par)
       if (allocated(par%x_ref)) deallocate(par%x_ref)

       if (par%multiprecision_used_internally) then 
#ifdef GMPXX_FOUND
         call free_pade_model(par%params_ptr)
#endif
       end if 

   end subroutine free_params

end module gx_ac
