! **************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
module test_pade_approximant
  ! External libs
  use zofu, only: unit_test_type
  ! Our libs/modules
  use kinds, only: dp

  ! Module being tested
  use gx_ac, only: thiele_pade_api, create_thiele_pade, evaluate_thiele_pade_at, & 
                   free_params, params
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


  !> @brief helper function to assert if the difference between two complex 
  !!        numbers is under a threshold
  !!
  !! Unfortunately Zofu asserts require real or double precision declarations
  !! (etc), which are not equivalent to real(kind=dp), complex(kind=dp)...
  !! so one needs a helper function such that a logical can be evaluated.
  !! Sensible thing would be to fork the framework and modify, or open a PR.
  !!
  !! @param[in] a   - complex number
  !! @param[in] b   - complex number
  !! @param[in] tol - threshold
  !! @return (abs(a - b) <= tol)
  logical function is_close_complex_dp(a, b, tol)
    complex(kind=dp), intent(in) :: a
    complex(kind=dp), intent(in) :: b
    real(kind=dp), optional, intent(in) :: tol
    ! abs() evaluates to real, hence tolerance is real
    real(kind=dp) :: tolerance = 1.-8_dp
    if (present(tol)) tolerance = tol
    is_close_complex_dp = abs(a - b) <= tolerance
  end function is_close_complex_dp



  !> @brief Test the Pade interpolant against the function -1 / (x - x0)
  !!
  !! @param test - test object
  subroutine test_pade(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Variable and function, respectively
    complex(kind=dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and the its reference value
    complex(kind=dp) :: f_approx, ref
    !> Some function center
    complex(kind=dp), parameter :: x0 = cmplx(2.0_dp, 2.0_dp, kind=dp)
    complex(kind=dp), parameter :: xx = cmplx(1.0_dp, 1.0_dp, kind=dp)
    !> Tolerance
    real(kind=dp) :: tol = 1.e-7_dp
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



  !> @brief Test the GMP Pade interpolant against the function -1 / (x - x0)
  !!
  !! @param test - test object
  subroutine test_pade_mp(test)
    !> Test object
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    type(params) :: params_thiele
    !> Variable and function, respectively
    complex(kind=dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and its reference value
    complex(kind=dp), dimension(:), allocatable :: xx
    complex(kind=dp), dimension(:), allocatable :: f_approx
    complex(kind=dp) :: ref
    !> Some function center
    complex(kind=dp), parameter :: x0 = cmplx(2.0_dp, 2.0_dp, kind=dp)
    !> Tolerance
    real(kind=dp) :: tol = 1.e-7_dp
    integer :: i

    !> Test setup
    allocate(x(n), f(n), xx(1), f_approx(1))
    do i = 1, n
       x(i) = cmplx(i, 0.1_dp, kind=dp)
       f(i) = -1.0_dp / (x(i) - x0)
    end do

    xx = cmplx(1.5_dp, 0.1_dp, kind=dp)
    ref = -1.0_dp / (xx(1) - x0)

    params_thiele = create_thiele_pade(n, x, f)
    f_approx(1:1) = evaluate_thiele_pade_at(params_thiele, xx)

    !> Test execution
    call test%assert(is_close(f_approx(1), ref, tol=tol), name = 'Test Pade GMP ~ -1 / (x - x0)')

    !> Clean-up
    deallocate(x, f, xx, f_approx)
    call free_params(params_thiele)

  end subroutine test_pade_mp



  !> @brief Test the Thiele-Pade interpolant against the function 1 / (-x^2 + 1) 
  !!        which has poles
  !!
  !! @param test - test object
  subroutine test_thiele_pade_poles(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Variable, function, and parameters, respectively
    complex(kind=dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and its reference value
    complex(kind=dp) :: ref
    complex(kind=dp), dimension(1) :: f_approx
    !> Test point
    complex(kind=dp), dimension(1), parameter :: xx = cmplx(1.0_dp, 3.0_dp, kind=dp)
    !> Tolerance
    real(kind=dp) :: tol = 1.e-7_dp
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



  !> @brief Test the GMP Thiele-Pade interpolant against the function 
  !!        1 / (-x^2 + 1) which has poles
  !!
  !! @param test - test object
  subroutine test_thiele_pade_poles_mp(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Variable, function, and parameters, respectively
    complex(kind=dp), allocatable :: x(:), f(:)
    type(params)          :: params_thiele
    !> Pade approximant of f, and its reference value
    complex(kind=dp), dimension(:), allocatable :: xx 
    complex(kind=dp), dimension(:), allocatable :: f_approx
    complex(kind=dp) :: ref
    !> Tolerance
    real(kind=dp) :: tol = 1.e-7_dp
    integer :: i

    !> Test setup
    allocate(x(n), f(n), xx(1), f_approx(1))
    do i = 1, n
       x(i) = cmplx(i, 0.05_dp, kind=dp)
       f(i) = 1.0_dp / (-x(i) * x(i) + 1.0_dp)
    end do

    xx = cmplx(1.0_dp, 3.0_dp, kind=dp) 
    ref = 1.0_dp / (-xx(1) * xx(1) + 1.0_dp)
    params_thiele = create_thiele_pade(n, x, f)
    f_approx(1:1) = evaluate_thiele_pade_at(params_thiele, xx)

    !> Test execution
    call test%assert(is_close(f_approx(1), ref, tol=tol), name = 'Test GMP Thiele-Pade ~ 1 / (-x^2 + 1)')

    !> Clean-up
    deallocate(x, f, xx, f_approx)
    call free_params(params_thiele)

  end subroutine test_thiele_pade_poles_mp



  !> @brief Test the Thiele-Pade interpolant against the function |x| which has 
  !!        a branch point
  !!
  !! @param test - test object
  subroutine test_thiele_pade_abs(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Newman grid constant
    real(kind=dp), parameter :: eta = exp(-1.0_dp / sqrt(dble(n)))
    real(kind=dp), parameter :: delta_eta = 0.0005_dp
    !> Variable, function, and parameters, respectively
    complex(kind=dp), allocatable :: x(:), f(:)
    !> Pade approximant of f, and its reference value
    complex(kind=dp) :: ref
    complex(kind=dp), dimension(1) :: f_approx
    !> Test point
    complex(kind=dp), dimension(1), parameter :: xx = cmplx(0.7_dp, 0.0_dp, kind=dp)
    !> Tolerance
    real(kind=dp) :: tol = 1.e-7_dp
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



  !> @brief Test the GMP Thiele-Pade interpolant against the function |x| which 
  !!        has a branch point
  !!
  !! @param test - test object
  subroutine test_thiele_pade_abs_mp(test)
    class(unit_test_type), intent(inout) :: test

    !> N sampling points
    integer, parameter :: n = 100
    !> Newman grid constant
    real(kind=dp), parameter :: eta = exp(-1.0_dp / sqrt(dble(n)))
    real(kind=dp), parameter :: delta_eta = 0.0005_dp
    !> Variable, function, and parameters, respectively
    complex(kind=dp), allocatable :: x(:), f(:)
    type(params)          :: params_thiele
    !> Pade approximant of f, and its reference value
    complex(kind=dp), dimension(:), allocatable :: xx 
    complex(kind=dp), dimension(:), allocatable :: f_approx
    complex(kind=dp) :: ref
    !> Test point
    !> Tolerance
    real(kind=dp) :: tol = 1.e-7_dp
    integer :: i, npar

    !> Test setup
    npar = 2 * n
    allocate(x(npar), f(npar), xx(1), f_approx(1))

    !> Here we use a Newman grid with 2n points
    do i = 1, n
       x(i) = cmplx(-eta**(i - 1), 0.0_dp, kind=dp)
       x(i+n) = cmplx(eta**(n - i) + delta_eta, 0.0_dp, kind=dp)
    end do

    f(:) = abs(x(:))

    xx = cmplx(0.7_dp, 0.0_dp, kind=dp)
    ref = abs(xx(1))

    params_thiele = create_thiele_pade(npar, x, f)
    f_approx(1:1) =  evaluate_thiele_pade_at(params_thiele, xx)

    !> Test execution
    call test%assert(is_close(f_approx(1), ref, tol=tol), name = 'Test GMP Thiele-Pade ~ |x|')

    !> Clean-up
    deallocate(x, f, xx, f_approx)
    call free_params(params_thiele)

  end subroutine test_thiele_pade_abs_mp

end module test_pade_approximant
