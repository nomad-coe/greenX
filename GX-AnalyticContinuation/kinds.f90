!> Precision constants
module kinds
   implicit none

   public
   integer, parameter :: dp=kind(1.0d0)                ! double presicion for real number
   integer, parameter :: dpc=kind((1.0_dp,1.0_dp))     ! double presicion for complex number
   real(dp), parameter :: zero=0._dp                   ! real type zero 
   complex(dpc), parameter :: czero = (0._dp,0._dp)    ! complex type zero
   complex(dpc), parameter :: cone  = (1._dp,0._dp)    ! complex type one

end module kinds
