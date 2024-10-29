!> @brief natural cubic splines functionality
!!
!! implementation according to 
!! Numerical Recipes, 3rd edition, cambridge university press, page 120
module spline 

    use kinds, only: dp 

    implicit none 

    private 
    public :: cubic_spline

    !> @brief natural cubic spline to interpolate 1D functions
    type cubic_spline 
        !> number of reference points
        integer :: n 
        !> grid of radii
        real(kind=dp), dimension(:), allocatable :: r_grid 
        !> tabulated function values
        real(kind=dp), dimension(:), allocatable :: func 
        !> tabulated second derivative
        real(kind=dp), dimension(:), allocatable :: sec_dev

        contains 

        procedure :: create => create_cubic_spline
        procedure :: evaluate => evaluate_cubic_spline
        procedure :: evaluate_batch => evaluate_cubic_spline_batch
        procedure :: finalize => deallocate_cubic_spline_arrays

    end type cubic_spline
    

    contains 

    !> @brief set up the natural cubic spline 
    !!
    !!        for natural splines the second derivative at the grid boundaries 
    !!        is zero 
    !!
    !! @param[in] r_grid -- grid of radius 
    !! @param[in] func   -- function values 
    subroutine create_cubic_spline(this, r_grid, func) 
        class(cubic_spline), intent(inout) :: this
        real(kind=dp), dimension(:) :: r_grid
        real(kind=dp), dimension(:) :: func 

        this%n = size(r_grid)
        if (this%n .ne. size(func)) then 
            print *, "ERROR: create_cubic_spline: r_grid and func not the same length!"
            stop 
        end if 
        allocate(this%r_grid(this%n), this%func(this%n), this%sec_dev(this%n))
        this%r_grid(:) = r_grid(:)
        this%func(:) = func(:)

        call calc_second_derivative(this)

    end subroutine create_cubic_spline

    

    !> @brief evaluate the spline at given r 
    !!
    !! @param[in] r -- radius
    real(kind=dp) function evaluate_cubic_spline(this, r) result(y)
        class(cubic_spline), intent(in) :: this 
        real(kind=dp), intent(in) :: r 

        ! internal variables
        integer :: klo, khi 
        real(kind=dp) :: h, b, a 

        klo = get_x_idx_in_array(r, this%r_grid)
        if (klo .eq. 0) then 
            print *, "ERROR: evaluate_cubic_spline: r query either outside of splined function or equal to spline grid"
        end if 
        khi = klo + 1
        h = this%r_grid(khi) - this%r_grid(klo)
        a = (this%r_grid(khi) - r)/h 
        b = (r-this%r_grid(klo))/h 
        y = a * this%func(klo) + b* this%func(khi) &
            + ((a**3 - a) * this%sec_dev(klo) + (b**3 - b) * this%sec_dev(khi)) &
            * h**2 / 6.0_dp

    end function evaluate_cubic_spline



    !> @brief evaluate the spline for batch of radii
    !!
    !! @param[in] r -- radii
    !! @param[out] y -- evaluated function values 
    subroutine evaluate_cubic_spline_batch(this, r, y_out)
        class(cubic_spline), intent(in) :: this 
        real(kind=dp), dimension(:), intent(in) :: r
        real(kind=dp), dimension(:), intent(out) :: y_out

        ! internal variables
        integer :: n, i

        n = size(r)
        do i = 1, n 
            y_out(i) = this%evaluate(r(i))
        end do 

    end subroutine evaluate_cubic_spline_batch



    !> @brief deallocate arrays of cubic_spline struct
    subroutine deallocate_cubic_spline_arrays(this)
        class(cubic_spline), intent(inout) :: this 
        deallocate(this%r_grid)
        deallocate(this%func)
        deallocate(this%sec_dev)
    end subroutine
    
    

    !> @brief get the second derivative if the function
    subroutine calc_second_derivative(splined)
        type(cubic_spline), intent(inout) :: splined

        ! internal variables 
        integer :: i, k 
        real(kind=dp) :: p, qn, sig, un 
        real(kind=dp), dimension(splined%n-1) :: u

        splined%sec_dev(1) = 0.0_dp
        u(1) = 0.0_dp
        
        do i = 2, splined%n -1
            sig = (splined%r_grid(i) - splined%r_grid(i-1)) &
                  / (splined%r_grid(i+1)-splined%r_grid(i-1))
            p = sig * splined%sec_dev(i-1) + 2.0_dp 
            splined%sec_dev(i) = (sig - 1.0_dp)/p 
            u(i) = (splined%func(i+1) - splined%func(i)) &
                   / (splined%r_grid(i+1) - splined%r_grid(i)) &
                   - (splined%func(i) - splined%func(i-1)) &
                   / (splined%r_grid(i) - splined%r_grid(i-1))
            u(i) = (6.0_dp * u(i)/(splined%r_grid(i+1) - splined%r_grid(i-1)) &
                   - sig*u(i-1))/p
        end do 

        qn = 0.0_dp
        un = 0.0_dp

        splined%sec_dev(splined%n) = (un - qn * u(splined%n-1)) &
                                / (qn * splined%sec_dev(splined%n-1) + 1.0_dp)
        do k = splined%n-1, 1, -1
            splined%sec_dev(k) = splined%sec_dev(k)*splined%sec_dev(k+1) + u(k)
        end do 

    end subroutine calc_second_derivative



    !> @brief find the first index where (x > x_array(idx)) and (x < x_array(idx+1))
    !!
    !! a return of 0 signals that x is outside of the array values or is 
    !! exactly equal to one element of the array
    !!
    !! @param[in] x       -- radius
    !! @param[in] x_array -- array of radii
    !! @return               index of the array
    integer function get_x_idx_in_array(x, x_array) result(idx)
        real(kind=dp), intent(in) :: x 
        real(kind=dp), dimension(:), intent(in) :: x_array 

        ! internal variables 
        integer :: n, i

        n = size(x_array)
        if ((x .le. x_array(1)) .or. (x > x_array(n))) then 
            idx = 0
            return 
        end if 
        
        n = size(x_array)
        do i = 1, n - 1
            if (x .eq. x_array(i + 1)) then 
                idx = 0
                return 
            else if  (x < x_array(i + 1)) then 
                idx = i
                return 
            end if 
        end do 
        idx = 0

    end function get_x_idx_in_array

end module spline 
