!> @brief math helper functions
module math

    use kinds, only: dp

    implicit none 

    private 
    public :: cart_to_sphere, double_factorial

    contains

    !> @brief convert a cartesian coordinate to a spherical coordinate 
    !!
    !! @param[in] cart -- cartesian coordinate
    !! @param[out] r -- radius
    !! @param[out] theta -- polar angle 
    !! @param[out] phi -- azimuth
    subroutine cart_to_sphere(cart, r, theta, phi)
        real(kind=dp), dimension(3), intent(in) :: cart
        real(kind=dp), intent(out) :: r       
        real(kind=dp), intent(out) :: theta  
        real(kind=dp), intent(out) :: phi
        
        ! Compute the radial distance
        r = sqrt(cart(1)**2 + cart(2)**2 + cart(3)**2)
        ! Compute the polar angle (theta)
        if (r > 0.0d0) then
            theta = acos(cart(3) / r)
        else
            theta = 0.0d0  ! Define theta as 0 if the vector has zero length
        end if
        ! Compute the azimuthal angle (phi)
        if (cart(1) == 0.0d0 .and. cart(2) == 0.0d0) then
            phi = 0.0d0  ! Define phi as 0 if x and y are both zero
        else
            phi = atan2(cart(2), cart(1))
            !phi = sign(1.0d0, y)*acos(x/sqrt(x**2 + y**2))
        end if
    end subroutine cart_to_sphere



    !> @brief returns the double fatcorial n!! of an integer number
    !!
    !! @param[in] n -- positive integer
    !! @return         double factorial 
    real(kind=dp) function double_factorial(n) result(dfact)
        integer, intent(in) :: n

        ! internal variables
        integer :: k 

        if (n < 0) then 
            print *, "ERROR: double factorial of negative number"
            stop 
        end if 

        dfact = 1.0_dp 
        if (n.eq.1 .or. n.eq.0) return

        if (mod(n, 2)==0) then 
            do k = 1, n/2
                dfact = dfact * (2.0_dp*k)
            end do 
        else 
            do k = 1, (n+1)/2
                dfact = dfact * (2.0_dp*k - 1.0_dp)
            end do
        end if 

    end function double_factorial

end module math 