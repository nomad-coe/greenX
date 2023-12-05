! **************************************************************************************************
!  Copyright (C) 2020-2023 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
module test_pade_approximant
  ! External libs
  use zofu, only: unit_test_type
  ! Our libs/modules
  use kinds, only: dp

  ! Module being tested
  use gx_ac, only: thiele_pade_api, create_thiele_pade_mp, evaluate_thiele_pade_mp, params_mp
  use pade_approximant, only: pade

  implicit none
  private
  public :: test_pade, test_thiele_pade_poles, test_thiele_pade_abs, &
            test_pade_mp, test_thiele_pade_poles_mp, test_thiele_pade_abs_mp

  ! Helper function
  interface is_close
     module procedure is_close_complex_dp
  end interface is_close

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
  end function is_close_complex_dp

  !> Test the Pade interpolant against the function -1 / (x - x0)
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
    complex(dp), parameter :: x0 = cmplx(2.0_dp, 2.0_dp, kind=dp)
    complex(dp), parameter :: xx = cmplx(1.0_dp, 1.0_dp, kind=dp)
    !> Tolerance
    real(dp) :: tol = 1.e-7_dp
    integer :: i

    !> Test setup
    allocate(x(n), f(n))
    do i = 1, n
       x(i) = cmplx(i, 0.0_dp, kind=dp)
       f(i) = -1.0_dp / (x(i) - x0)
    end do

    ref = cmplx(0.5, -0.5, dp)
    f_approx = pade(n, x, f, xx)

    !> Test execution
    call test%assert(is_close(f_approx, ref, tol=tol), name = 'Test Pade ~ -1 / (x - x0)')

    !> Clean-up
    deallocate(x)
    deallocate(f)

  end subroutine test_pade

  !> Test the GMP Pade interpolant against the function -1 / (x - x0)
  subroutine test_pade_mp(test)
    !> Test object
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    type(params_mp) :: params
    !> Variable and function, respectively
    complex(dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and its reference value
    complex(dp) :: xx
    complex(dp) :: f_approx
    complex(dp) :: ref
    !> Some function center
    complex(dp), parameter :: x0 = cmplx(2.0_dp, 2.0_dp, kind=dp)
    !> Tolerance
    real(dp) :: tol = 1.e-7_dp
    integer :: i

    !> Test setup
    allocate(x(n), f(n))
    do i = 1, n
       x(i) = cmplx(i, 0.0_dp, kind=dp)
       f(i) = -1.0_dp / (x(i) - x0)
    end do

    xx = cmplx(1.0_dp, 1.0_dp, kind=dp)
    ref = -1.0_dp / (xx - x0)

    call create_thiele_pade_mp(n, x, f, params)
    call evaluate_thiele_pade_mp(xx, f_approx, params)


    !> Test execution
    call test%assert(is_close(f_approx, ref, tol=tol), name = 'Test Pade GMP ~ -1 / (x - x0)')

    !> Clean-up
    deallocate(x)
    deallocate(f)

  end subroutine test_pade_mp

  !> Test the Thiele-Pade interpolant against the function 1 / (-x^2 + 1) which has poles
  subroutine test_thiele_pade_poles(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Variable, function, and parameters, respectively
    complex(dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and its reference value
    complex(dp) :: ref
    complex(dp), dimension(1) :: f_approx
    !> Test point
    complex(dp), dimension(1), parameter :: xx = cmplx(1.0_dp, 3.0_dp, kind=dp)
    !> Tolerance
    real(dp) :: tol = 1.e-7_dp
    integer :: i

    !> Test setup
    allocate(x(n), f(n))
    do i = 1, n
       x(i) = cmplx(i, 0.05_dp, kind=dp)
       f(i) = 1.0_dp / (-x(i) * x(i) + 1.0_dp)
    end do
    ref = 1.0_dp / (-xx(1) * xx(1) + 1.0_dp)

    call thiele_pade_api(n, x, f, xx, f_approx, .true.)

    !> Test execution
    call test%assert(is_close(f_approx(1), ref, tol=tol), name = 'Test Thiele-Pade ~ 1 / (-x^2 + 1)')

    !> Clean-up
    deallocate(x)
    deallocate(f)

  end subroutine test_thiele_pade_poles

  !> Test the GMP Thiele-Pade interpolant against the function 1 / (-x^2 + 1) which has poles
  subroutine test_thiele_pade_poles_mp(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Variable, function, and parameters, respectively
    complex(dp), allocatable :: x(:), f(:)
    type(params_mp)          :: params
    !> Pade approximant of f, and its reference value
    complex(dp):: xx 
    complex(dp) :: f_approx
    complex(dp) :: ref
    !> Tolerance
    real(dp) :: tol = 1.e-7_dp
    integer :: i

    !> Test setup
    allocate(x(n), f(n))
    do i = 1, n
       x(i) = cmplx(i, 0.05_dp, kind=dp)
       f(i) = 1.0_dp / (-x(i) * x(i) + 1.0_dp)
    end do

    xx = cmplx(1.0_dp, 3.0_dp, kind=dp) 
    ref = 1.0_dp / (-xx * xx + 1.0_dp)
    call create_thiele_pade_mp(n, x, f, params)
    call evaluate_thiele_pade_mp(xx, f_approx, params)

    !> Test execution
    call test%assert(is_close(f_approx, ref, tol=tol), name = 'Test GMP Thiele-Pade ~ 1 / (-x^2 + 1)')

    !> Clean-up
    deallocate(x)
    deallocate(f)

  end subroutine test_thiele_pade_poles_mp

  !> Test the Thiele-Pade interpolant against the function |x| which has a branch point
  subroutine test_thiele_pade_abs(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Newman grid constant
    real(dp), parameter :: eta = exp(-1.0_dp / sqrt(dble(n)))
    real(dp), parameter :: delta_eta = 0.0005_dp
    !> Variable, function, and parameters, respectively
    complex(dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and its reference value
    complex(dp) :: ref
    complex(dp), dimension(1) :: f_approx
    !> Test point
    complex(dp), dimension(1), parameter :: xx = cmplx(0.7_dp, 0.0_dp, kind=dp)
    !> Tolerance
    real(dp) :: tol = 1.e-7_dp
    integer :: i, npar

    !> Test setup
    npar = 2 * n
    allocate(x(npar), f(npar))

    !> Here we use a Newman grid with 2n points
    do i = 1, n
       x(i) = cmplx(-eta**(i - 1), 0.0_dp, kind=dp)
       x(n + i) = cmplx(eta**(n - i) + delta_eta, 0.0_dp, kind=dp)
    end do

    f(:) = abs(x(:))
    ref = abs(xx(1))

    call thiele_pade_api(npar, x, f, xx, f_approx, .true.)

    !> Test execution
    call test%assert(is_close(f_approx(1), ref, tol=tol), name = 'Test Thiele-Pade ~ |x|')

    !> Clean-up
    deallocate(x)
    deallocate(f)

  end subroutine test_thiele_pade_abs

  !> Test the GMP Thiele-Pade interpolant against the function |x| which has a branch point
  subroutine test_thiele_pade_abs_mp(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Newman grid constant
    real(dp), parameter :: eta = exp(-1.0_dp / sqrt(dble(n)))
    real(dp), parameter :: delta_eta = 0.0005_dp
    !> Variable, function, and parameters, respectively
    complex(dp), allocatable :: x(:), f(:)
    type(params_mp)          :: params
    !> Pade approximant of f, and its reference value
    complex(dp) :: xx 
    complex(dp) :: f_approx
    complex(dp) :: ref
    !> Test point
    !> Tolerance
    real(dp) :: tol = 1.e-7_dp
    integer :: i, npar

    !> Test setup
    npar = 2 * n
    allocate(x(npar), f(npar))

    !> Here we use a Newman grid with 2n points
    do i = 1, n
       x(i) = cmplx(-eta**(i - 1), 0.0_dp, kind=dp)
       x(i+n) = cmplx(eta**(n - i) + delta_eta, 0.0_dp, kind=dp)
    end do

    f(:) = abs(x(:))

    xx = cmplx(0.7_dp, 0.0_dp, kind=dp)
    ref = abs(xx)

    call create_thiele_pade_mp(npar, x, f, params)
    call evaluate_thiele_pade_mp(xx, f_approx, params)


    !> Test execution
    call test%assert(is_close(f_approx, ref, tol=tol), name = 'Test GMP Thiele-Pade ~ |x|')

    !> Clean-up
    deallocate(x)
    deallocate(f)

  end subroutine test_thiele_pade_abs_mp

end module test_pade_approximant
