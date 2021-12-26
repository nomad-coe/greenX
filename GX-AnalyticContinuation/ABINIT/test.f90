! COPYRIGHT
!  Copyright (C) 2020-2021 Green-X library (MA, MG, XG)
!  This file is distributed under the terms of the
!  APACHE2 License.

module test

  use kinds
  use pade_approximant

  implicit none

  public :: test1
  public :: test2

contains

! FUNCTION test1
! Produce the results of pade approximant for a given function, in this case f(x) = 1 / w - w0, 
! while w and w0 are real and imaginary values respectively. 

function test1()

    complex(dpc):: test1
    integer, parameter :: n = 100
    complex(dpc), parameter :: w = (2.0,2.0)
    complex(dpc),parameter :: xx = (1.0,1.0) 

    complex(dpc), parameter :: j = (1._dp,1._dp)
  
    complex(dpc) :: f(n)
    complex(dpc) :: x(n)
    integer :: i, m
    i = 1
    do i = 1,n
        x(i) = i
        f(i) = -1/(x(i)-w) 
       
    enddo
    test1 = pade(n,x,f,xx)
 
end function 
 
!----------------------------------------------------------------------

! FUNCTION test2
! Produce the results of pade approximant for a given function, in this case f(x) = 1 / w - w0, 
! while w and w0 are real and imaginary values respectively. 

function test2()

    complex(dpc):: test2
    integer, parameter :: n = 100
    complex(dpc), parameter :: w = (2.0,2.0)
    complex(dpc),parameter :: xx = (1.0,1.0) 

    complex(dpc), parameter :: j = (1._dp,1._dp)
  
    complex(dpc) :: f(n)
    complex(dpc) :: x(n)
    integer :: i, m
    i = 1
    do i = 1,n
        x(i) = i
        f(i) = -1/sqrt(x(i)-w) 
       
    enddo
    test2 = pade(n,x,f,xx)
 
end function 
 
!----------------------------------------------------------------------

end module test
