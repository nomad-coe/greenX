! ***************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

!>  Fortran double precision implementation of the thiele pade approximation. 
module pade_approximant
   use kinds, only: dp
   implicit none

   private
   public :: pade, pade_derivative, thiele_pade, evaluate_thiele_pade, evaluate_thiele_pade_tab

   !> Complex zero
   complex(dp), public :: c_zero = cmplx(0.0_dp, 0.0_dp, kind=dp)
   !> Complex one
   complex(dp), public :: c_one = cmplx(1.0_dp, 0.0_dp, kind=dp)

contains

   !> @brief Calculate the pade approximant in $xx$ point of the function $f_n(x)$
   !> calculated at the $n$ points $x$
   !>
   !> @param[in]  n  Number of points
   !> @param[in]  x  Variable evaluated at discrete points {n}
   !> @param[in]  f  Function to approximate
   !> @param[in]  xx  Pade will be computed for this value
   !> @return     pade   Pade approximant
   complex(dp) function pade(n, x, f, xx)
      integer, intent(in)     :: n
      complex(dp), intent(in) :: xx
      complex(dp), intent(in) :: x(n), f(n)

      !> Pade coefficients
      complex(dp) :: a(n)

      ! Generate parameters
      call pade_coefficient_derivative(x, f, a)

      ! Evaluate using Wallis method
      call evaluate_thiele_pade(n, x, xx, a, pade, "none")

   end function pade

   !> @brief Calculate the derivative of the pade approximant in xx of the
   !> function f calculated at the n points x.
   !>
   !> @param[in]  n  Number of points
   !> @param[in]  x  Variable evaluated at discrete points {n}
   !> @param[in]  f  Function to approximate
   !> @param[in]  xx  Pade will be computed for this value
   !> @return     pade   Derivative of the pade approximant
   complex(dp) function pade_derivative(n, x, f, xx)
      integer, intent(in)     :: n
      complex(dp), intent(in) :: xx
      complex(dp), intent(in) :: x(n), f(n)

      integer :: i
      complex(dp) :: a(n)
      !> Coefficients in the numerator and denominator, respectively
      complex(dp) :: acoef(0:n), bcoef(0:n)
      !> Derivatives are acoef and bcoef
      complex(dp) :: dacoef(0:n), dbcoef(0:n)

      call pade_coefficient_derivative(x, f, a)

      acoef(0) = c_zero
      acoef(1) = a(1)
      bcoef(0:1) = c_one
      dacoef(0:1) = c_zero
      dbcoef(0:1) = c_zero

      do i = 1, n - 1
         acoef(i + 1) = acoef(i) + (xx - x(i))*a(i + 1)*acoef(i - 1)
         bcoef(i + 1) = bcoef(i) + (xx - x(i))*a(i + 1)*bcoef(i - 1)
         dacoef(i + 1) = dacoef(i) + a(i + 1)*acoef(i - 1) + (xx - x(i))*a(i + 1)*dacoef(i - 1)
         dbcoef(i + 1) = dbcoef(i) + a(i + 1)*bcoef(i - 1) + (xx - x(i))*a(i + 1)*dbcoef(i - 1)
      end do
      pade_derivative = dacoef(n)/bcoef(n) - acoef(n)*dbcoef(n)/(bcoef(n)*bcoef(n))

   end function pade_derivative

   !> @brief Calculate the derivative of the the coefficients of pade approximant
   !>
   !> @param[in]  x  Variable evaluated at discrete points {n}
   !> @param[in]  f  Function to approximate
   !> @return     a   Derivative of the pade approximant coefficients
   subroutine pade_coefficient_derivative(x, f, a)
      complex(dp), intent(in)  :: x(:), f(:)
      complex(dp), intent(out) :: a(:)

      ! Internal variables
      integer                  :: i, j, n
      complex(dp), allocatable :: c(:, :)

      n = size(x)
      allocate (c(n, n))
      c(1, :) = f(:)

      do i = 2, n
         do j = i, n
            c(i, j) = (c(i - 1, i - 1) - c(i - 1, j))/((x(j) - x(i - 1))*c(i - 1, j))
         end do
      end do

      do i = 1, n
         a(i) = c(i, i)
      end do

   end subroutine pade_coefficient_derivative

   !> brief Gets the Pade approximant of a meromorphic function F
   !>       This routine implements a modified version of the Thiele's reciprocal differences
   !>       interpolation algorithm using a greedy strategy, ensuring that the ordering of the
   !>       included points minimizes the value of |P_n(x_{1+1}) - F(x_{i+1})|
   !>       The default Thiele interpolation is also included for conveniency
   !!  @param[in]  n_par - order of the interpolant
   !!  @param[inout] x_ref - array of the reference points
   !!  @param[in]  y_ref - array of the reference function values
   !!  @param[in]  do_greedy - whether to use the greedy algorithm or the naive one
   !!  @param[in]  enforce_symmetry - enforce given function symmetry (e.g. even)
   !!  @param[out] par - array of the interpolant parameters
   subroutine thiele_pade(n_par, x_ref, y_ref, a_par, do_greedy, enforce_symmetry)
      integer, intent(in)                           :: n_par
      complex(kind=dp), dimension(:), intent(inout) :: x_ref
      complex(kind=dp), dimension(:), intent(in)    :: y_ref
      complex(kind=dp), dimension(:), intent(out)   :: a_par
      character(*),  intent(in)                     :: enforce_symmetry
      logical, intent(in)                           :: do_greedy


      ! Internal variables
      integer                                       :: i, i_par, idx, jdx, kdx, n_rem
      integer, dimension(n_par)                     :: n_rem_idx
      real(kind=dp), parameter                      :: tol = 1.0E-6_dp
      real(kind=dp)                                 :: deltap, pval
      complex(kind=dp)                              :: pval_in, x_in, y_in, acoef_in, bcoef_in
      complex(kind=dp), dimension(n_par, n_par)     :: g_func
      complex(kind=dp), dimension(n_par)            :: x, xtmp, y_ref_tmp, ytmp
      complex(kind=dp), dimension(-1:n_par)         :: acoef, bcoef

      ! for some symmetry contraints y_ref changes in projection 
      do i_par = 1, n_par 
          y_ref_tmp(i_par) = apply_sym_y(x_ref(i_par), y_ref(i_par), enforce_symmetry)
      end do 

      ! Initialize variables
      acoef_in = c_zero
      bcoef_in = c_zero
      x_in = c_zero
      y_in = c_zero

      n_rem_idx(:) = (/(i, i = 1, n_par)/)
      a_par(:) = c_zero
      g_func(:, :) = c_zero
      x(:) = c_zero

      ! Whether to perform the refined Thiele's interpolation (default)
      if (do_greedy) then
         ! Unpack initial reference arguments, as they will be overwritten
         x(:) = x_ref
         x_ref = c_zero

         ! Select first point that maximizes |F|
         kdx = maxloc(abs(y_ref_tmp), dim=1)
         xtmp(1) = x(kdx)
         ytmp(1) = y_ref_tmp(kdx)
         x_ref(1) = x(kdx)

         n_rem = n_par - 1
         do i = kdx, n_rem
            n_rem_idx(i) = n_rem_idx(i + 1)
         end do

         ! Compute the generating function for the first time
         call thiele_pade_gcoeff(xtmp, ytmp, g_func, 1, enforce_symmetry)
         a_par(1) = g_func(1, 1)

         ! Initialize Walli's coefficients
         acoef(-1) = c_one
         acoef(0) = c_zero
         bcoef(-1) = c_zero
         bcoef(0) = c_one

         ! Add remaining points ensuring min |P_i(x_{1+1}) - F(x_{i+1})|
         do idx = 2, n_par
            pval = huge(0.0_dp)
            do jdx = 1, n_rem
               ! Compute next convergent P_i(x_{i+1})
               call evaluate_thiele_pade_tab(idx - 1, xtmp(1:idx-1), x(n_rem_idx(jdx)), &
                                             a_par, acoef, bcoef, enforce_symmetry)
               pval_in = acoef(idx - 1) / bcoef(idx - 1)

               ! Select the point that minimizes difference's absolute value
               deltap = abs(pval_in - y_ref_tmp(n_rem_idx(jdx)))
               if (deltap .lt. pval) then
                  pval = deltap
                  x_in = x(n_rem_idx(jdx))
                  y_in = y_ref_tmp(n_rem_idx(jdx))
                  acoef_in = acoef(idx - 1)
                  bcoef_in = bcoef(idx - 1)
                  kdx = jdx
               end if
            end do

            ! Update indexes of non-visited points
            n_rem = n_rem - 1
            do i = kdx, n_rem
               n_rem_idx(i) = n_rem_idx(i + 1)
            end do

            ! Add the winning point
            x_ref(idx) = x_in
            xtmp(idx) = x_in
            ytmp(idx) = y_in

            ! Rescale Wallis coefficients to avoid overflow
            acoef(idx - 1) = acoef_in
            bcoef(idx - 1) = bcoef_in
            if (abs(bcoef_in) > tol) then
               acoef(idx - 1) = acoef(idx - 1) / bcoef_in
               acoef(idx - 2) = acoef(idx - 2) / bcoef_in
               bcoef(idx - 1) = c_one
               bcoef(idx - 2) = bcoef(idx - 2) / bcoef_in
            end if

            ! Get the recurrence matrix
            call thiele_pade_gcoeff(xtmp, ytmp, g_func, idx, enforce_symmetry)

            ! Unpack parameters a_i = g_i(w_i)
            a_par(idx) = g_func(idx, idx)
         end do
      else
         ! Directly interpolate
         do i_par = 1, n_par
            call thiele_pade_gcoeff(x_ref, y_ref_tmp, g_func, i_par, enforce_symmetry)
            a_par(i_par) = g_func(i_par, i_par)
         enddo
      end if

   end subroutine thiele_pade

   !> brief Computes the recurrence coefficients from Thiele's continued fraction
   !>       This routine uses tabulation in order to efficienly compute the matrix elements g_func(:,:)
   !! @param[in] n - number of parameters
   !! @param[in] x - array of the reference points
   !! @param[in] y - array of the reference function values
   !! @param[in] symm - enforce function symmetry
   !! @param[inout] g_func - recurrence matrix used to compute the parameters a_n
   subroutine thiele_pade_gcoeff(x, y, g_func, n, symm)
      integer, intent(in)                             :: n
      complex(kind=dp), dimension(:), intent(in)      :: x, y
      complex(kind=dp), dimension(:,:), intent(inout) :: g_func
      character(len=*), intent(in)                    :: symm

      ! Internal variables
      integer :: idx

      ! Begin work (leveraging tabulation of the g_func)
      g_func(n, 1) = y(n)
      if (n==1) return

      do idx = 2, n
         g_func(n, idx) = (g_func(idx - 1, idx - 1) - g_func(n, idx - 1)) / &
            ((apply_sym_x(x(n), symm) - apply_sym_x(x(idx - 1), symm)) * g_func(n, idx - 1))
      enddo

   end subroutine thiele_pade_gcoeff

   !> brief Evaluates a Pade approximant constructed with Thiele's reciprocal differences
   !>       This is the tabulated version of the procedure
   !! @param[in] n_par - number of parameters
   !! @param[in] x_ref - array of the reference points
   !! @param[in] x - the point to evaluate
   !! @param[in] a_par -  array of the input parameters
   !! @param[in] acoef - parameter wallis algorithm
   !! @param[in] bcoef - parameter wallis algorithm
   !! @param[in] enforce_symm - enforce given symmetry
   !! @param[out] y -  the value of the interpolant at x
   subroutine  evaluate_thiele_pade_tab(n_par, x_ref , x, a_par, acoef, bcoef, symm)
      integer, intent(in)                             :: n_par
      complex(kind=dp), dimension(:), intent(in)      :: x_ref
      complex(kind=dp), intent(in)                    :: x
      complex(kind=dp), dimension(:), intent(in)      :: a_par
      complex(dp), dimension(-1:n_par), intent(inout) :: acoef, bcoef
      character(*), intent(in)                        :: symm

      ! Internal variables
      complex(dp)                                     :: delta

      ! Wallis' method iteration
      delta = a_par(n_par)
      if (n_par > 1) then
            delta = delta * (apply_sym_x(x, symm) - apply_sym_x(x_ref(n_par - 1), symm))
      end if

      acoef(n_par) = acoef(n_par - 1) + delta * acoef(n_par - 2)
      bcoef(n_par) = bcoef(n_par - 1) + delta * bcoef(n_par - 2)

   end subroutine evaluate_thiele_pade_tab

   !> brief Evaluates a Pade approximant constructed with Thiele's reciprocal differences
   !! @param[in] n_par - number of parameters
   !! @param[in] x_ref - array of the reference points
   !! @param[in] x - the point to evaluate
   !! @param[in] a_par -  array of the input parameters
   !! @param[out] y -  the value of the interpolant at x
   !! @param[in] enforce_symmetry - symmetry of the underlying function
   subroutine  evaluate_thiele_pade(n_par, x_ref , x, a_par, y, enforce_symmetry)
      integer, intent(in)                        :: n_par
      complex(kind=dp), dimension(:), intent(in) :: x_ref
      complex(kind=dp), intent(in)               :: x
      complex(kind=dp), dimension(:), intent(in) :: a_par
      complex(kind=dp), intent(out)              :: y
      character(*), intent(in)                   :: enforce_symmetry

      ! Internal variables
      integer                                    :: i_par
      real(kind=dp), parameter                   :: tol = 1.0E-6_dp
      complex(dp), dimension(-1:n_par)           :: acoef, bcoef

      ! Evaluate using Wallis' method
      acoef(-1) = c_one
      acoef(0) = c_zero
      bcoef(-1) = c_zero
      bcoef(0) = c_one

      do i_par = 1, n_par
         call evaluate_thiele_pade_tab(i_par, x_ref , x, a_par, acoef, &
                                       bcoef, enforce_symmetry)
         if (abs(bcoef(i_par)) > tol) then
            acoef(i_par) = acoef(i_par) / bcoef(i_par)
            acoef(i_par - 1) = acoef(i_par - 1) / bcoef(i_par)
            bcoef(i_par - 1) = bcoef(i_par - 1) / bcoef(i_par)
            bcoef(i_par) = c_one
         end if
      end do

      y = acoef(n_par) / bcoef(n_par)
      y = apply_sym_y(x, y, enforce_symmetry)

   end subroutine evaluate_thiele_pade

   !> brief for a given complex argument return the argument with enforced symmetry
   !!
   !! @param[in] x     - complex number 
   !! @param[in] symm  - symmetry
   !! @return    x_symm - symmetry constrained x
   complex(kind=dp) function apply_sym_x(x, symm) result(x_symm)
      complex(kind=dp), intent(in) :: x 
      character(len=*), intent(in) :: symm

      select case (symm)
        case ("y")
            x_symm = cmplx(abs(x%re), x%im, kind=8)
        case ("x")
            x_symm = cmplx(x%re, abs(x%im), kind=8)
        case ("xy")
            x_symm = cmplx(abs(x%re), abs(x%im), kind=8)
        case ("even")
            x_symm = cmplx(abs(x%re), sign(1.0d0, x%re) * x%im, kind=8)
        case ("odd")
            x_symm = cmplx(abs(x%re), sign(1.0d0, x%re) * x%im, kind=8)
        case ("conjugate")
            x_symm = cmplx(abs(x%re), sign(1.0d0, x%re) * x%im, kind=8)
        case ("anti-conjugate")
            x_symm = cmplx(abs(x%re), sign(1.0d0, x%re) * x%im, kind=8)
        case ("none")
            x_symm = x
        case default 
            print *, "*** thiele_pade: enforce_symmetry=", symm, &
                    " not known or not supported!" 
            stop 
      end select 

   end function apply_sym_x

   !> brief projects function argument y if x gets projected 
   !!
   !! @param[in] x_original - the function argument before projection
   !! @param[in] y - the function value before projection
   !! @param[in] symm  - symmetry constaint
   !! @return    y_symm - symmetry constrained y
   complex(kind=dp) function apply_sym_y(x_original, y, symm) result(y_symm)
      complex(kind=dp), intent(in) :: x_original
      complex(kind=dp), intent(in) :: y
      character(len=*), intent(in) :: symm

      ! internal variables
      complex(kind=dp):: x_projected

      x_projected = apply_sym_x(x_original, symm)

      if (x_projected.eq.x_original) then 
          ! no projection happened
          y_symm = y 
      else
          ! x got projected 
          select case (symm)
              case ("odd")
                  y_symm = -y
              case ("conjugate")
                  y_symm = conjg(y)
              case ("anti-conjugate")
                  y_symm = -conjg(y)
              case default 
                  y_symm = y 
          end select 
      end if 

   end function apply_sym_y

end module pade_approximant
