module log_grid 

    use kinds, only: dp 

    implicit none

    private 
    public :: create_log_grid

    contains 

    !> @brief create a logarithmically spaced grid
    !!
    !! the grid is constructed that the point before the first point is a r=0
    !! and the point behind the last point is at r=\infty 
    !!
    !! @param[in] r_outer -- outermost grid point 
    !! @param[in] N       -- number of grid points 
    !! @param[out] r_grid -- logarithmic radial grid
    subroutine create_log_grid(r_outer, N, r_grid)
        real(kind=dp), intent(in) :: r_outer 
        integer,       intent(in) :: N 
        real(kind=dp), dimension(N), intent(out) :: r_grid

        ! internal variables
        integer :: i 
        
        do i = 1, N 
            r_grid(i) = r_outer * log(1.0_dp - (i/(N+1.0_dp))**2)/log(1.0_dp - (N/(N+1.0_dp))**2)
        end do 

    end subroutine create_log_grid

end module log_grid