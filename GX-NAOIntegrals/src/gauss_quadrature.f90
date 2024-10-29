!> @brief gauss-quadrature using gauss-legendre grids
!!
!! see Numerical Recipes, 3rd edition, cambridge university press, page 184
module gauss_quadrature 

    use kinds, only: dp
    use constants, only: pi

    implicit none

    private 
    public :: get_gauss_legendre_grid, gauss_legendre_integrator

    !> precision for getting the gauss-legendre grid
    real(kind=dp), parameter :: gauleg_precision = 1.0e-14_dp

    contains 

    !> @brief get gauss-legendre abscissas and weights for integration
    !!
    !! see Numerical Recipes, 3rd edition, cambridge university press, page 184 
    !!
    !! @param[in] N         -- number of weights and abscissas
    !! @param[in] x_start   -- lower integral bound
    !! @param[in] x_stop    -- upper integral bound
    !! @param[out] x        -- tabulated abscissas (integration points)
    !! @param[out] w        -- tabulated weights
    subroutine get_gauss_legendre_grid(n, x_start, x_stop, x, w)
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: x_start
        real(kind=dp), intent(in) :: x_stop
        real(kind=dp), dimension(n), intent(out) :: x
        real(kind=dp), dimension(n), intent(out) :: w

        ! internal variables 
        real(kind=dp) :: z1, z, xm, xl, pp, p3, p2, p1
        integer :: m, i, j

        if (mod(n, 2) .ne. 0) then 
            print *, "ERROR: get_gauss_legendre_grid: odd value for n"
            stop 
        end if 

        m = (n+1)/2
        xm = 0.5_dp * (x_stop + x_start)
        xl = 0.5_dp * (x_stop - x_start)
        do i = 0, m - 1 
            z = cos(pi * (i + 0.75_dp)/(n+0.5_dp))
            z1 = huge(z)
            do while (abs(z - z1) > gauleg_precision)
                p1 = 1.0_dp 
                p2 = 0.0_dp 
                do j = 0, n-1
                    p3 = p2 
                    p2 = p1 
                    p1 = ((2.0_dp * j + 1.0_dp) * z * p2 - j*p3)/(j + 1.0_dp)
                end do 
                pp = n*(z*p1-p2) / (z**2 - 1.0_dp)
                z1 = z 
                z = z1 - p1/pp 
            end do 
            x(i+1) = xm - xl*z 
            x(n-i) = xm + xl*z 
            w(i+1) = 2.0_dp * xl / ((1.0_dp-z**2) * pp**2)
            w(n-i) = w(i+1)
        end do 

    end subroutine

    !> @brief evaluate gauss-legendre integral 
    !!
    !! @param[in] n -- integration grid points 
    !! @param[in] w -- gauss-legendre weights 
    !! @param[in] f -- function values on gauss-legendre grid 
    !! @return         value of the integral 
    real(kind=dp) function gauss_legendre_integrator(n, w, f) result(integral)
        integer, intent(in) :: n 
        real(kind=dp), dimension(n), intent(in) :: w
        real(kind=dp), dimension(n), intent(in) :: f

        integral = sum(w * f)

    end function gauss_legendre_integrator 


end module gauss_quadrature