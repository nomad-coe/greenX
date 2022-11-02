!  Copyright (C) 2020-2022 Green-X library
!  This file is distributed under the terms of the APACHE2 License.

! TODO(Alex) Move this documentation to below
!   The Pade approximants are a particular type of rational fraction
!   approximation to the value of a function. The idea is to match the Taylor
!   series expansion as far as possible.
!   Here, we Implemented the Pad\'e approximant using Thiele's reciprocal-difference method.
!   This routine takes a function $(f_n=f(x_n))$, considering complex $x_n$ which is
!   evaluated at an initial set of arguments, $(x_n)$
!   approximates the function with the help of Pad\'e approximants, and evaluates (extrapolates/rotates)
!   this approximation at a given set of arguments $(x)$. The $N$-point Pad\'e approximant
!   then reads
!   $$ P_N(x)=
!   \cfrac{a_1}
!   {1+\cfrac{a_2(x-x_1)}{\cdots+\cfrac{a_n(x-x_{N-1})}{1+(x-x_N)g_{N+1}(x)}}}
!   $$
!   \cfrac{a_1}
!   {1+\cfrac{a_2(x-x_1)}{\cdots+\cfrac{a_n(x-x_{N-1})}{1+(x-x_N)g_{N+1}(x)}}}
!   $$
!   where
!   $$  P_N(x)=
!   \lim_{n \to \infty}\cfrac{A_n(x)}{B_n(x)},
!   $$
!   $$  g_n(x)=\frac{g_{n-1}(x_{n-1})-g_{n-1}(x)}
!                   {(x-x_{n-1})g_{n-1}(x)}, \; n \ge 2
!   $$
!   and
!   $$  a_n=g_n(x_n),\; g_1(x_n)=f_n,\; n=1,\ldots,N.
!   $$
!
!   Expressions are taken from G. A. J. Baker, Essentials of Padé Approximants (Academic,New York, 1975).
!   See also PHYSICAL REVIEW B 94, 165109 (2016).
module pade_approximant
   use kinds, only: dp
   implicit none

   private
   public :: pade, pade_derivative

   !> Complex zero
   complex(dp) :: c_zero = cmplx(0._dp, 0._dp, kind=dp)
   !> Complex one
   complex(dp) :: c_one = cmplx(1._dp, 0._dp, kind=dp)

contains

   ! TODO(Alex) Replace this with description from above
   !> @brief Calculate the pade approximant in $xx$ point of the function $f_n(x)$
   !> calculated at the $n$ points $x$
   !>
   !> @param[in]  n  Number of points
   !> @param[in]  x  ADD ME
   !> @param[in]  f  ADD ME
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
         acoef(i + 1) = acoef(i) + (xx - x(i))*a(i + 1)*acoef(i - 1)
         bcoef(i + 1) = bcoef(i) + (xx - x(i))*a(i + 1)*bcoef(i - 1)
      end do

      pade = acoef(n)/bcoef(n)

   end function pade

   !> @brief Calculate the derivative of the pade approximant in xx of the
   !> function f calculated at the n points x.
   !>
   !> @param[in]  n  Number of points
   !> @param[in]  x   ADD ME
   !> @param[in]  f   ADD ME
   !> @param[in]  xx  ADD ME
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
   !> @param[in]  x   ADD ME
   !> @param[in]  f   ADD ME
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

end module pade_approximant
