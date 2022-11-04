module test_pade_approximant
    ! External libs
    use zofu, only: unit_test_type
    ! Our libs/modules
    use kinds, only: dp

    ! Module being tested
    use pade_approximant, only: pade

    implicit none
    private
    public :: test_pade

    ! Helper function
    interface is_close
        module procedure is_close_complex_dp
    end interface

contains

    ! Unfortunately Zofu asserts require real or double precision declarations
    ! (etc), which are not equivalent to real(dp), complex(dp)...
    ! so one needs a helper function such that a logical can be evaluated.
    ! Sensible thing would be to fork the framework and modify, or open a PR.
    logical function is_close_complex_dp(a, b, tol)
        complex(dp), intent(in) :: a
        complex(dp), intent(in) :: b
        real(dp), optional, intent(in) :: tol
        ! abs() evaluates to real, hence tolerance is real
        real(dp) :: tolerance = 1.-8_dp
        if (present(tol)) tolerance = tol
        is_close_complex_dp = abs(a - b) <= tolerance
    end function

    !> Test pade against the function -1 / (x - x0)
    subroutine test_pade(test)
        !> Test object
        class(unit_test_type), intent(inout) :: test

        !> N sampling points
        integer, parameter :: n = 100
        !> Variable and function, respectively
        complex(dp), allocatable :: x(:), f(:)
        !> Pade approximant of f, and the its reference value
        complex(dp) :: f_approx, ref
        !> Some function center
        complex(dp), parameter :: x0 = (2.0, 2.0)
        complex(dp), parameter :: xx = (1.0, 1.0)
        !> Tolerance
        real(dp) :: tol = 1.e-7_dp
        integer :: i

        allocate(x(n), f(n))
        do i = 1, n
            x(i) = cmplx(i, 0, kind=dp)
            f(i) = -1._dp / (x(i) - x0)
        end do

        ref = cmplx(0.5, -0.5, dp)
        f_approx = pade(n, x, f, xx)

        call test%assert(is_close(f_approx, ref, tol=tol), name = 'Test pade ~ -1 / (x - x0)')

    end subroutine test_pade

end module test_pade_approximant
