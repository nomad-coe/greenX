! ***************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

!>   The Pade approximants are a particular type of rational fraction
!>   approximation to the value of a function. The idea is to match the Taylor
!>   series expansion as far as possible.
!>   Here, we Implemented the Pad\'e approximant using Thiele's reciprocal-difference method.
!>   This routine takes a function $(f_n=f(x_n))$, considering complex $x_n$ which is
!>   evaluated at an initial set of arguments, $(x_n)$
!>   approximates the function with the help of Pad\'e approximants, and evaluates (extrapolates/rotates)
!>   this approximation at a given set of arguments $(x)$. The $N$-point Pad\'e approximant
!>   then reads
!>   $$ P_N(x)=
!>     \cfrac{a_1}
!>     {1+\cfrac{a_2(x-x_1)}{\cdots+\cfrac{a_n(x-x_{N-1})}{1+(x-x_N)g_{N+1}(x)}}}
!>   $$
!>     \cfrac{a_1}
!>     {1+\cfrac{a_2(x-x_1)}{\cdots+\cfrac{a_n(x-x_{N-1})}{1+(x-x_N)g_{N+1}(x)}}}
!>   $$
!>   where
!>   $$  P_N(x)=
!>          \lim_{n \to \infty}\cfrac{A_n(x)}{B_n(x)},
!>   $$
!>   $$  g_n(x)=\frac{g_{n-1}(x_{n-1})-g_{n-1}(x)}
!>                   {(x-x_{n-1})g_{n-1}(x)}, \; n \ge 2
!>   $$
!>   and
!>   $$  a_n=g_n(x_n),\; g_1(x_n)=f_n,\; n=1,\ldots,N.
!>   $$
!>
!>   Expressions are taken from G. A. J. Baker, Essentials of Padé Approximants (Academic,New York, 1975).
!>   See also:
!>   PHYSICAL REVIEW B 94, 165109 (2016).
!>   J. CHEM. THEORY COMPUT. 19, 16, 5450–5464 (2023)
module pade_approximant
   use kinds, only: dp
   implicit none

   private
   public :: pade, pade_derivative, thiele_pade, evaluate_thiele_pade

   !> Complex zero
   complex(dp) :: c_zero = cmplx(0._dp, 0._dp, kind=dp)
   !> Complex one
   complex(dp) :: c_one = cmplx(1._dp, 0._dp, kind=dp)

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
      integer, intent(in) :: n
      complex(dp), intent(in) :: xx
      complex(dp), intent(in) :: x(n), f(n)

      !> Pade coefficients
      complex(dp) :: a(n)
      !> Numerator and denominator in pade series
      complex(dp) :: acoef(0:n), bcoef(0:n)
      integer :: i

      call pade_coefficient_derivative(x, f, a)

      acoef(0) = c_zero
      acoef(1) = a(1)
      bcoef(0:1) = c_one

      do i = 1, n - 1
         acoef(i + 1) = acoef(i) + (xx - x(i)) * a(i + 1) * acoef(i - 1)
         bcoef(i + 1) = bcoef(i) + (xx - x(i)) * a(i + 1) * bcoef(i - 1)
      end do

      pade = acoef(n)/bcoef(n)

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
      integer, intent(in) :: n
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
      complex(dp), intent(in) :: x(:), f(:)
      complex(dp), intent(out) :: a(:)

      integer :: i, j, n
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
   !!  @param[in]  do_greedy - whether to use the default greedy algorithm or the naive one
   !!  @param[out] par - array of the interpolant parameters
   subroutine thiele_pade(n_par, x_ref, y_ref, a_par, do_greedy)
      integer, intent(in)                            :: n_par
      complex(kind=dp), dimension(:), intent(inout)  :: x_ref
      complex(kind=dp), dimension(:), intent(in)     :: y_ref
      complex(kind=dp), dimension(:), intent(inout)  :: a_par
      logical, optional, intent(in)                  :: do_greedy

      ! Internal variables
      logical                                        :: local_do_greedy = .True.
      integer                                        :: i, i_par, idx, jdx, kdx, n_rem
      integer, dimension(n_par)                      :: n_rem_idx
      real(kind=dp)                                  :: deltap, pval
      complex(kind=dp)                               :: pval_in, x_in, y_in
      complex(kind=dp), dimension(n_par, n_par)      :: g_func
      complex(kind=dp), dimension(n_par)             :: x, xtmp, ytmp

      ! Whether to perform the refined Thiele's interpolation (default)
      if (present(do_greedy)) local_do_greedy = do_greedy

      ! Initialize arrays
      n_rem_idx = (/(i, i = 1, n_par)/)
      g_func = c_zero
      x_in = c_zero
      y_in = c_zero

      if (local_do_greedy) then
         ! Unpack initial reference arguments, as they will be overwritten
         x = x_ref
         x_ref = c_zero

         ! Select first point as to maximize |Wmn|
         kdx = minloc(abs(x),dim=1)
         xtmp(1) = x(kdx)
         ytmp(1) = y_ref(kdx)
         x_ref(1) = x(kdx)

         n_rem = n_par - 1
         do i = kdx, n_rem
            n_rem_idx(i) = n_rem_idx(i + 1)
         end do

         ! Compute the generating function for the first time
         call thiele_pade_gcoeff(xtmp, ytmp, g_func, 1)
         a_par(1) = g_func(1, 1)

         ! Add remaining points ensuring min |P_i(x_{1+1}) - Wmn(x_{i+1})|
         do idx = 2, n_par
            pval = huge(0.0_dp)
            do jdx = 1, n_rem
               ! Compute next convergent P_i(x_{i+1})
               call evaluate_thiele_pade(idx - 1, xtmp(1:idx-1), x(n_rem_idx(jdx)), a_par, pval_in)

               ! Select the point that minimizes difference's absolute value
               deltap = abs(pval_in - y_ref(n_rem_idx(jdx)))
               if (deltap .lt. pval) then
                  pval = deltap
                  x_in = x(n_rem_idx(jdx))
                  y_in = y_ref(n_rem_idx(jdx))
                  kdx = jdx
               end if
            end do

            ! Update indexes of non-visited points
            n_rem = n_rem - 1
            do i = kdx, n_rem
               n_rem_idx(i) = n_rem_idx(i+1)
            end do

            ! Add the winning point and recompute generating function
            x_ref(idx) = x_in
            xtmp(idx) = x_in
            ytmp(idx) = y_in
            call thiele_pade_gcoeff(xtmp,ytmp,g_func,idx)

            ! Unpack parameters a_i = g_i(w_i)
            do i_par = 1, idx
               a_par(i_par) = g_func(i_par, i_par)
            enddo
         end do
      else
         ! Interpolate
         call thiele_pade_gcoeff(x_ref, y_ref, g_func, 1)
         a_par(1) = g_func(1, 1)
         do i_par = 2, n_par
            call thiele_pade_gcoeff(x_ref, y_ref, g_func, i_par)
            a_par(i_par) = g_func(i_par, i_par)
         enddo
      end if

   end subroutine thiele_pade

   !> brief Computes the recurrence coefficients from Thiele's continued fraction
   !>       This routine uses tabulation in order to efficienly compute the matrix elements g_func(:,:)
   !! @param[in] n - number of parameters
   !! @param[in] x - array of the frequencies
   !! @param[in] y - array of the Wmn matrix elements
   !! @param[inout] g_func - recurrence matrix used to compute the parameters a_n
   subroutine thiele_pade_gcoeff(x, y, g_func, n)
      integer, intent(in)                             :: n
      complex(kind=dp), dimension(:), intent(in)      :: x, y
      complex(kind=dp), dimension(:,:), intent(inout) :: g_func

      ! Internal variables
      integer :: idx

      ! Begin work (leveraging tabulation of the g_func)
      g_func(n, 1) = y(n)
      if (n==1) return

      do idx = 2, n
         g_func(n, idx) = (g_func(idx - 1, idx - 1) - g_func(n, idx - 1)) / &
            ((x(n) - x(idx - 1)) * (x(n) + x(idx - 1)) * g_func(n, idx - 1))
      enddo
   end subroutine thiele_pade_gcoeff

   !> brief Gets the value of the Wmn matrices using the previously computed Pade approximant
   !>       Here we only implement the Pade approximant evaluation
   !! @param[in] n_par - number of parameters
   !! @param[in] x_ref - array of the reference points
   !! @param[in] x - the point to evaluate
   !! @param[in] a_par -  array of the input parameters
   !! @param[out] y -  the value of the interpolant at x
   subroutine  evaluate_thiele_pade(n_par, x_ref , x, a_par, y)
      integer, intent(in)                        :: n_par
      complex(kind=dp), dimension(:), intent(in) :: x_ref
      complex(kind=dp), intent(in)               :: x
      complex(kind=dp), dimension(:), intent(in) :: a_par
      complex(kind=dp), intent(inout)            :: y

      ! Internal variables
      integer                                    :: i_par
      complex(kind=dp)                           :: gtmp

      ! Begin work
      gtmp = c_one

      do i_par = n_par, 2, -1
         gtmp = 1.0_dp + a_par(i_par) * (x - x_ref(i_par - 1)) * (x + x_ref(i_par - 1)) / gtmp
      enddo

      ! Compute out value
      y = a_par(1) / gtmp

   end subroutine evaluate_thiele_pade

end module pade_approximant
