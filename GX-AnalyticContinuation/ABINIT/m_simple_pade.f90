module m_simple_pade

  implicit none

  public :: pade

  public
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
  real(dp), parameter :: zero=0._dp
  complex(dpc), parameter :: czero = (0._dp,0._dp)
  complex(dpc), parameter :: cone  = (1._dp,0._dp)



contains

!----------------------------------------------------------------------

!!****f* m_numeric_tools/pade
!! NAME
!!  pade
!!
!! FUNCTION
!!  Calculate the pade approximant in zz of the function f calculated at the n points z
!!
!! SOURCE

function pade(n,z,f,zz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: zz
 complex(dpc) :: pade
!arrays
 complex(dpc),intent(in) :: z(n),f(n)

!Local variables-------------------------------
!scalars
 complex(dpc) :: a(n)
 complex(dpc) :: Az(0:n), Bz(0:n)
 integer :: i
! *************************************************************************

 call calculate_pade_a(a,n,z,f)

 Az(0)=czero
 Az(1)=a(1)
 Bz(0)=cone
 Bz(1)=cone

 do i=1,n-1
   Az(i+1)=Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
   Bz(i+1)=Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
   !write(*,*) 'A & B =', Az(i+1), Bz(i+1)
 end do
 !write(std_out,*) 'Bz(n)',Bz(n)
 !if (REAL(Bz(n))==zero.and.AIMAG(Bz(n))==zero) write(std_out,*) ' Bz(n) ',Bz(n)
 pade=Az(n)/Bz(n)
 !write(*,*) 'pade =', pade
 !write(std_out,*) 'pade_approx ', pade_approx

end function pade
!!***
function test()
    complex(dpc):: test
    integer, parameter :: n = 300
    complex(dpc),parameter :: zz = (0.1,0.0) 

    complex(dpc), parameter :: j = (0._dp,1._dp)
  
    complex(dpc) :: f(n)
    complex(dpc) :: z(n)
    integer :: i, m
    i = 1
    do i = 1,n


        z(i) = i *j

        f(i) = cos(z(i)) 
        
        
        
        !write(*,*) f(i)
    enddo
  
    !test = exp(zz)
    test = pade(n,z,f,zz)
    !write(*,*) test
  
end function 
!----------------------------------------------------------------------

!!****f* m_numeric_tools/dpade
!! NAME
!!  dpade
!!
!! FUNCTION
!!  Calculate the derivative of the pade approximant in zz of the function f calculated at the n points z
!!
!! SOURCE

function dpade(n,z,f,zz)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: zz
 complex(dpc) :: dpade
!arrays
 complex(dpc),intent(in) :: z(n),f(n)

!Local variables-------------------------------
!scalars
 integer :: i
!arrays
 complex(dpc) :: a(n)
 complex(dpc) :: Az(0:n), Bz(0:n)
 complex(dpc) :: dAz(0:n), dBz(0:n)
! *************************************************************************

 call calculate_pade_a(a,n,z,f)

 Az(0)=czero
 Az(1)=a(1)
 Bz(0)=cone
 Bz(1)=cone
 dAz(0)=czero
 dAz(1)=czero
 dBz(0)=czero
 dBz(1)=czero

 do i=1,n-1
   Az(i+1)=Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
   Bz(i+1)=Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
   dAz(i+1)=dAz(i)+a(i+1)*Az(i-1)+(zz-z(i))*a(i+1)*dAz(i-1)
   dBz(i+1)=dBz(i)+a(i+1)*Bz(i-1)+(zz-z(i))*a(i+1)*dBz(i-1)
 end do
 !write(std_out,*) 'Bz(n)', Bz(n)
 !if (REAL(Bz(n))==zero.and.AIMAG(Bz(n))==zero) write(std_out,*) 'Bz(n)',Bz(n)
 !pade_approx = Az(n) / Bz(n)
 dpade=dAz(n)/Bz(n) -Az(n)*dBz(n)/(Bz(n)*Bz(n))
 !write(std_out,*) 'pade_approx ', pade_approx

end function dpade
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/calculate_pade_a
!! NAME
!!  calculate_pade_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_numeric_tools
!!
!! CHILDREN
!!
!! SOURCE

subroutine calculate_pade_a(a,n,z,f)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 complex(dpc),intent(in) :: z(n),f(n)
 complex(dpc),intent(out) :: a(n)

!Local variables-------------------------------
!scalars
 integer :: i,j
!arrays
 complex(dpc) :: g(n,n)
! *************************************************************************

 g(1,1:n)=f(1:n)

 do i=2,n
   do j=i,n
     !if (REAL(g(i-1,j))==zero.and.AIMAG(g(i-1,j))==zero) write(std_out,*) 'g_i(z_j)',i,j,g(i,j)
     g(i,j)=(g(i-1,i-1)-g(i-1,j)) / ((z(j)-z(i-1))*g(i-1,j))
     !write(std_out,*) 'g_i(z_j)',i,j,g(i,j)
   end do
 end do
 do i=1,n
   a(i)=g(i,i)
 end do
 !write(std_out,*) 'a ',a(:)

end subroutine calculate_pade_a
!!***





end module m_simple_pade
