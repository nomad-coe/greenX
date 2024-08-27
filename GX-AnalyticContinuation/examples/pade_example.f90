! ***************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

!> @brief example of using the AC-component to do analytic continuation 
!!        of a two-pole model function
!!
!! This example program is intendet to showcase the usage of the GreenX - Analytic 
!! Continuation component by fitting a 16-parameter pade model at a two-pole function.
!!
!!      reference function:     y = (x - 0.25)^{-1} + (x - 0.75)^{-1} - 50
!!
!! - the reference function values are taken along the imaginary axis x \in [0i, 1i]
!! - these values are used to obtain a pade model 
!! - this pade model is used to evaluate the function along the real axis with 
!!   a small imaginary shift x \in [0 + imag_shift, 1 + imag_shift]
!! - the function values obtained with the pade model are compared against the 
!!   correct ones
program use_pade 

    use gx_ac, only: create_thiele_pade, evaluate_thiele_pade_at, & 
                     free_params, params

    implicit none 

    integer,      parameter :: n_parameters = 16            ! number of pade parameters (= points on imaginary axis)
    integer,      parameter :: n_points     = 100           ! number of interpolated points on real axis
    real(kind=8), parameter :: imag_shift   = 0.001d0       ! imaginary shift for interpolated points 

    type(params)                               :: params_thiele
    complex(kind=8), dimension(:), allocatable :: x_query
    complex(kind=8), dimension(:), allocatable :: y_return
    complex(kind=8), dimension(:), allocatable :: y_correct
    complex(kind=8), dimension(:), allocatable :: x_ref
    complex(kind=8), dimension(:), allocatable :: y_ref

    allocate(x_ref(n_parameters), y_ref(n_parameters))
    allocate(x_query(n_points), y_return(n_points), y_correct(n_points)) 

    ! create points along the imaginary axis
    call create_complex_grid(n_parameters, x_ref, constant_along="real", shift=0.0d0)
    ! evaluate function values along imaginary axis 
    call evaluate_two_pole_model(n_parameters, x_ref, y_ref)
    ! create a grid along the real axis where analytically continued function 
    !    values will be evaluated by pade
    call create_complex_grid(n_points, x_query, constant_along="imag", shift=imag_shift)

    ! ---- GreenX API calls ------------------------------------------------------------

    ! create the pade interpolation model and store it in struct
    params_thiele = create_thiele_pade(n_parameters, x_ref, y_ref, &
                                       do_greedy=.false., precision=64)

    ! evaluate the pade interpolation model at given x points
    y_return(1:n_points) =  evaluate_thiele_pade_at(params_thiele, x_query)

    ! ----------------------------------------------------------------------------------

    ! compare the analytic continuation with the correct function values
    call evaluate_two_pole_model(n_points, x_query, y_correct)
    call print_comparison(n_points, x_query, y_return, y_correct)


    ! Clean-up
    call free_params(params_thiele)
    deallocate(x_query)
    deallocate(y_return)
    deallocate(y_correct)
    deallocate(x_ref)
    deallocate(y_ref)

    contains 
        ! definition of some helper functions

        !> @brief get function values of a two-pole model for a given complex 
        !!        grid
        !!
        !! @parameter[in]  n_points -- number of points
        !! @parameter[in]  x        -- coplex grid points
        !! @parameter[out] y        -- function values of two-pole model
        subroutine evaluate_two_pole_model(n_points, x, y)
            integer,                              intent(in)  :: n_points 
            complex(kind=8), dimension(n_points), intent(in)  :: x
            complex(kind=8), dimension(n_points), intent(out) :: y

            y = (x - 0.25d0)**(-1.0d0) + (x - 0.75d0)**(-1.0d0) - 50.0d0*x

        end subroutine evaluate_two_pole_model


        !> @brief creates a complex grid along imaginary or real axis
        !! 
        !! @parameter[in]  n_points       -- number of grid points
        !! @parameter[in]  constant_along -- axis that is constant 
        !! @parameter[out] x              -- grid points along a given axis
        subroutine create_complex_grid(n_points, x, constant_along, shift)
            integer,                              intent(in)  :: n_points 
            complex(kind=8), dimension(n_points), intent(out) :: x
            character(*),                         intent(in)  :: constant_along
            real(kind=8),                         intent(in)  :: shift

            ! internal variables
            integer :: i

            real(kind=8), parameter :: x_start = 0.0d0
            real(kind=8), parameter :: x_end   = 1.0d0
            real(kind=8) :: x_step 
            real(kind=8), dimension(n_points) :: x_tmp

            complex(kind=8) :: prototype
            complex(kind=8) :: constant_prototype

            x_step = (x_end - x_start)/ real(n_points - 1, kind=8)
            do i = 1, n_points 
                x_tmp(i) = x_start + (i-1) * x_step 
            end do 

            if (constant_along .eq. "real") then 
                prototype = cmplx(0.0d0, 1.0d0, kind=8)
                constant_prototype = cmplx(shift, 0.0d0, kind=8)
            else if (constant_along .eq. "imag") then
                prototype = cmplx(1.0d0, 0.0d0, kind=8)
                constant_prototype = cmplx(0.0d0, shift, kind=8)
            end if 

            do i = 1, n_points 
                x(i) = x_tmp(i) * prototype + constant_prototype
            end do 

        end subroutine create_complex_grid


        !> @brief print complex (x,y) point pairs and for comparison the 
        !!        correct value y_ref
        !!
        !! @parameter[in] n_points -- number of points
        !! @parameter[in] x        -- complex function argument
        !! @parameter[in] y        -- complex function value interpolated by pade
        !! @parameter[in] y_ref    -- correct complex function value
        subroutine print_comparison(n_points, x, y_pade, y_ref)
            integer,                              intent(in)  :: n_points 
            complex(kind=8), dimension(n_points), intent(in)  :: x
            complex(kind=8), dimension(n_points), intent(out) :: y_pade
            complex(kind=8), dimension(n_points), intent(out) :: y_ref

            ! internal variables 
            integer :: i 

            print "(a, 2x, a, 2x, a, 2x, a, 2x, a, 2x, a)", "  real(x) ", & 
                  "imag(x)   ", "real(y_pade)  ", "imag(y_pade)  ", &
                  "real(y_ref)   ", "imag(y_ref)    "

            do i = 1, n_points
                print "(f10.4, 2x, f10.4, 2x, f14.8, 2x, f14.8, 2x, f14.8, 2x, f14.8)", &
                      x(i)%re, x(i)%im, y_pade(i)%re, y_pade(i)%im, y_ref(i)%re, y_ref(i)%im
            end do

        end subroutine print_comparison

end program use_pade