! ***************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

!>   The Pade approximants are a particular type of rational fraction
!!   approximation to the value of a function. The idea is to match the Taylor
!!   series expansion as far as possible.
!!   Here, we Implemented the Pade approximant using Thiele's reciprocal-difference method.
!!   This routine takes a function \f$f_n=f(x_n)\f$, considering complex \f$x_n\f$ which is
!!   evaluated at an initial set of arguments, \f$(x_n)\f$
!!   approximates the function with the help of Pade approximants, and evaluates (extrapolates/rotates)
!!   this approximation at a given set of arguments \f$(x)\f$. The \f$N\f$-point Pade approximant
!!   then reads
!!   $$ f(x) \approx P_N(x)=
!!     \cfrac{a_1}
!!     {1+\cfrac{a_2(x-x_1)}{\cdots+\cfrac{a_n(x-x_{N-1})}{1+(x-x_N)g_{N+1}(x)}}}
!!   $$
!!   where
!!   $$  g_n(x)=\frac{g_{n-1}(x_{n-1})-g_{n-1}(x)}
!!                   {(x-x_{n-1})g_{n-1}(x)}, \; n \ge 2
!!   $$
!!   and
!!   \f[  a_n=g_n(x_n)\\ g_1(x_n)=f_n\\ n=1,\ldots,N \f]
!!
!!   Expressions are taken from G. A. J. Baker, Essentials of Padé Approximants (Academic,New York, 1975).
!!   See also:
!!   PHYSICAL REVIEW B 94, 165109 (2016);
!!   J. CHEM. THEORY COMPUT. 19, 16, 5450–5464 (2023)
module gx_ac
   use kinds, only: dp
   use, intrinsic :: iso_c_binding, only: c_int, c_double_complex, c_ptr
   use pade_approximant, only: evaluate_thiele_pade, thiele_pade, c_zero, c_one
   implicit none

   private 
   public :: params, &
             create_thiele_pade, &
             evaluate_thiele_pade_at, &
             free_params, &
             arbitrary_precision_available, &
             thiele_pade_api


   !> @brief store the parameters of the thiele pade model 
   !!        (potentially in abitrary precision floats using GMP)
   type :: params 
       !> switch to check whether parameters are already initialized
       logical    :: initialized = .false.          
       !> number of pade parameters
       integer    :: n_par
       !> internal float arithmetic precision
       integer    :: precision 
       !> switch to check whether GMP was used for multi precision floats
       logical    :: multiprecision_used_internally
       !> enforce symmetry of pade fit e.g. ("x", "y", "xy", "even", "none")
       character(len=15) :: enforced_symmetry
       !> switch to check wether greedy algorithm was used 
       logical    :: use_greedy

       !> pointer to c++ struct for multiple precision arithmetic
       type(c_ptr) :: params_ptr

       !> reference points (used in fortran double precision routines)
       complex(kind=dp), dimension(:), allocatable :: x_ref
       !> pade parameters (used in fortran double precision routines)
       complex(kind=dp), dimension(:), allocatable :: a_par
   end type params 

#ifdef GMPXX_FOUND
   logical, parameter :: arbitrary_precision_available = .true.
#else
   logical, parameter :: arbitrary_precision_available = .false.
#endif

#ifdef GMPXX_FOUND
   interface

      !> @brief auxiliary function to compute Thiele-Pade parameters using 
      !!        arbitrary precision numbers
      !!
      !! @param[in] n_par - order of the interpolant
      !! @param[in] x_ref - array of the reference points
      !! @param[in] y_ref - array of the reference function values
      !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
      !! @return - pointer to abstract type to store all parameters
      function thiele_pade_mp_aux(n_par, x_ref, y_ref, do_greedy, precision, symmetry) bind(C, name="thiele_pade_mp")
         import :: c_double_complex, c_int, c_ptr
         integer(c_int), value                    :: n_par
         complex(c_double_complex), dimension(*)  :: x_ref
         complex(c_double_complex), dimension(*)  :: y_ref
         integer(c_int), value                    :: do_greedy
         integer(c_int), value                    :: precision
         integer(c_int), value                    :: symmetry
         type(c_ptr)                              :: thiele_pade_mp_aux
      end function thiele_pade_mp_aux

      !> @brief auxiliary function to evaluate the Thiele-Pade parameters using 
      !!        arbitrary precision numbers
      !!
      !! @param[in] x - point where the function needs to be evaluated
      !! @param[in] params_ptr - pointer to abstract type to store all parameters
      !! @return - interpolated function value
      function evaluate_thiele_pade_mp_aux(x, params_ptr) bind(C, name="evaluate_thiele_pade_mp")
         import :: c_ptr, c_double_complex
         complex(c_double_complex), value    :: x
         type(c_ptr), value                  :: params_ptr
         complex(c_double_complex)           :: evaluate_thiele_pade_mp_aux
      end function evaluate_thiele_pade_mp_aux

      !> @brief Frees the C++ pointers used in the Fortran type
      !!
      !! @param[in] params_ptr - the C++ pointer
      subroutine free_pade_model(params_ptr) bind(C, name="free_pade_model")
         import :: c_ptr
         type(c_ptr), value :: params_ptr
      end subroutine free_pade_model

   end interface
#endif

contains

   !> @brief API function to compute Thiele-Pade approximations of a meromorphic 
   !!        function
   !!       
   !!        >> Deprecated, please use create_thiele_pade() and 
   !!        evaluate_thiele_pade_at() whenever possible! <<
   !!
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
      call thiele_pade(n_par, x_ref_local, y_ref, a_par, do_greedy, "none")

      ! Compute the number of query points
      num_query = size(x_query)

      ! Evaluate the Thiele-Pade approximation at the query points
      do i = 1, num_query
         call evaluate_thiele_pade(n_par, x_ref_local, x_query(i), a_par, y_query(i), "none")
      end do

   end subroutine thiele_pade_api



   !> @brief API function to compute Thiele-Pade parameters 
   !!        (potentially using arbitrary precision arithmetic)
   !!
   !! @param[in] n_par - order of the interpolant
   !! @param[in] x     - array of the reference points
   !! @param[in] y     - array of the reference function values
   !! @param[in] do_greedy - whether to use the default greedy algorithm or the naive one
   !! @param[in] precision - precision in bits (!! not bytes !!)
   !! @param[in] enforce_symmetry - force the model to have a certain symmetry ("x", "y", "xy", "even", "none")
   !! @return    params - abstract type to store all parameters in arb. prec. representation
   type(params) function create_thiele_pade(n_par, x, y, do_greedy, precision, enforce_symmetry) result(par)
      integer, intent(in)                        :: n_par
      complex(kind=dp), dimension(:), intent(in) :: x, y
      logical, optional, intent(in)              :: do_greedy
      integer, optional, intent(in)              :: precision 
      character(*), optional, intent(in)         :: enforce_symmetry

      ! Internal variables
      integer :: c_do_greedy
      integer :: c_symmetry_label

      ! initialize type
      par%initialized = .true.
      par%n_par = n_par

      ! precision of internal arithmetic
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

      ! Symmetry consistency check
      if (present(enforce_symmetry)) then 
          if ((enforce_symmetry.eq."x") &
              .or. (enforce_symmetry.eq."y") & 
              .or. (enforce_symmetry.eq."xy") &
              .or. (enforce_symmetry.eq."even") &
              .or. (enforce_symmetry.eq."odd") &
              .or. (enforce_symmetry.eq."conjugate") &
              .or. (enforce_symmetry.eq."anti-conjugate") &
              .or. (enforce_symmetry.eq."none")) then
              ! character for fortran
              par%enforced_symmetry = enforce_symmetry
              ! integer for c
              select case (enforce_symmetry)
                  case ("y")
                      c_symmetry_label = 1
                  case ("x")
                      c_symmetry_label = 2
                  case ("xy")
                      c_symmetry_label = 3
                  case ("even")
                      c_symmetry_label = 4
                  case ("odd")
                      c_symmetry_label = 5
                  case ("conjugate")
                      c_symmetry_label = 6
                  case ("anti-conjugate")
                      c_symmetry_label = 7
                  case ("none")
                      c_symmetry_label = 0
                  case default 
                      print *, "symmetry not known!"
                      stop 
              end select 
          else
              print *, "*** create_thiele_pade: enorce_symmetry=", enforce_symmetry, &
                       " not known or not supported! Aborting..."
              stop
          end if  
      else 
          par%enforced_symmetry = "none"
          c_symmetry_label = 0
      end if 

      ! greedy algorithm
      if (present(do_greedy)) then 
          ! actual bool for fortran
          par%use_greedy = do_greedy
          ! use integer bools for interoperability with C
          if (do_greedy) then 
            c_do_greedy = 1
          else 
            c_do_greedy = 0
          end if 
      else 
          par%use_greedy = .true.
          c_do_greedy = 0
      end if 

      ! create the pade model 
      if (.not. par%multiprecision_used_internally) then 
          allocate(par%a_par(n_par))
          allocate(par%x_ref(size(x)))
          par%x_ref(:) = x
          call thiele_pade(n_par, par%x_ref, y, par%a_par, &
                           par%use_greedy, par%enforced_symmetry)
      elseif (par%multiprecision_used_internally) then 
#ifdef GMPXX_FOUND
          par%params_ptr = thiele_pade_mp_aux(n_par, x, y, c_do_greedy, &
                                              par%precision, c_symmetry_label)
#endif
      end if 

   end function create_thiele_pade 



   !> @brief API function to evaluate the Thiele-Pade approximation 
   !!        (potentially using arbitrary precision numbers)
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
            call evaluate_thiele_pade(par%n_par, par%x_ref, x(i), par%a_par, &
                                      y(i), par%enforced_symmetry)
        elseif (par%multiprecision_used_internally) then 
#ifdef GMPXX_FOUND
            y(i) = evaluate_thiele_pade_mp_aux(x(i), par%params_ptr)
#endif
        end if 

      end do

   end function evaluate_thiele_pade_at



   !> @brief deallocate the pade model 
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
