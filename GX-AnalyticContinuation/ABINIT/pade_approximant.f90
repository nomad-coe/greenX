! COPYRIGHT
!  Copyright (C) 2020-2021 Green-X library (MA, MG, XG)
!  This file is distributed under the terms of the
!  APACHE2 License.

module pade_approximant

  use kinds

! !DESCRIPTION:
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
! 

  implicit none

  public :: pade

contains

! FUNCTION test1
! Produce the results of pade approximant for a given function, in this case f(x) = 1 / w - w0, 
! while w and w0 are real and imaginary values respectively. 


!----------------------------------------------------------------------

! FUNCTION pade
!  Calculate the pade approximant in $xx$ point of the function $f_n(x)$ calculated at the $n$ points $x$
! 

function pade(n,x,f,xx)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n                            ! number of points
 complex(dpc),intent(in) :: xx
 complex(dpc) :: pade
!arrays
 complex(dpc),intent(in) :: x(n),f(n)               ! 

!Local variables-------------------------------
!scalars
 complex(dpc) :: a(n)                               ! pade coefficients
 complex(dpc) :: acoef(0:n), bcoef(0:n)                   ! numerator and denominator in pade series
 integer :: i

 call calculate_pade_coef(a,n,x,f)

 acoef(0)=czero                                        ! see kinds.f90
 acoef(1)=a(1)                                         ! initialization
 bcoef(0)=cone                                         ! see kinds.f90
 bcoef(1)=cone                                         ! see kinds.f90

 do i=1,n-1
   acoef(i+1)=acoef(i)+(xx-x(i))*a(i+1)*acoef(i-1)
   bcoef(i+1)=bcoef(i)+(xx-x(i))*a(i+1)*bcoef(i-1)
 end do
 pade=acoef(n)/bcoef(n)

end function pade

!----------------------------------------------------------------------

! FUNCTION dpade
!  Calculate the derivative of the pade approximant in xx of the function f calculated at the n points x


function dpade(n,x,f,xx)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n                            ! number of points
 complex(dpc),intent(in) :: xx                      ! pade will be calculated for this value
 complex(dpc) :: dpade
!arrays
 complex(dpc),intent(in) :: x(n),f(n)

!Local variables-------------------------------
!scalars
 integer :: i
!arrays
 complex(dpc) :: a(n)
 complex(dpc) :: acoef(0:n), bcoef(0:n)            ! acoef and bcoef are the coefficients in the nominator and denominator respectively 
 complex(dpc) :: dacoef(0:n), dbcoef(0:n)          ! derivatives are acoef and bcoef

 call calculate_pade_coef(a,n,x,f)

 acoef(0)=czero
 acoef(1)=a(1)
 bcoef(0)=cone
 bcoef(1)=cone
 dacoef(0)=czero
 dacoef(1)=czero
 dbcoef(0)=czero
 dbcoef(1)=czero

 do i=1,n-1
   acoef(i+1)=acoef(i)+(xx-x(i))*a(i+1)*acoef(i-1)
   bcoef(i+1)=bcoef(i)+(xx-x(i))*a(i+1)*bcoef(i-1)
   dacoef(i+1)=dacoef(i)+a(i+1)*acoef(i-1)+(xx-x(i))*a(i+1)*dacoef(i-1)
   dbcoef(i+1)=dbcoef(i)+a(i+1)*bcoef(i-1)+(xx-x(i))*a(i+1)*dbcoef(i-1)
 end do
 dpade=dacoef(n)/bcoef(n) -acoef(n)*dbcoef(n)/(bcoef(n)*bcoef(n))

end function dpade

!----------------------------------------------------------------------

! SUBROUTINE calculate_pade_coef
!  Calculate the derivative of the the coefficients of pade approximant 

subroutine calculate_pade_coef(a,n,x,f)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: x(n),f(n)
 complex(dpc),intent(out) :: a(n)

!Local variables-------------------------------
!scalars
 integer :: i,j
!arrays
 complex(dpc) :: c(n,n)

 c(1,1:n)=f(1:n)

 do i=2,n
   do j=i,n
     c(i,j)=(c(i-1,i-1)-c(i-1,j)) / ((x(j)-x(i-1))*c(i-1,j))
   end do
 end do
 do i=1,n
   a(i)=c(i,i)
 end do

end subroutine calculate_pade_coef

!----------------------------------------------------------------------

end module pade_approximant
