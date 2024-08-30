! ***************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************

!> @brief tests the thiele pade approximation using a model function 
!!          
!!      usage: 
!!                  ./test_gx_analytic_continuation <no_greedy/greedy:str> <precision:int>
!!      e.g. non-greedy algorithm and double precision: 
!!                  ./test_gx_analytic_continuation no_greedy 64
!!
!!      reference function:     y = (x - 0.25)^{-1} + (x - 0.75)^{-1} - 50
!!
!! - the reference function values are taken along the imaginary axis x \in [0i, 1i]
!! - these values are used to obtain a pade model 
!! - this pade model is used to evaluate the function along the real axis with 
!!   a small imaginary shift x \in [0 + imag_shift, 1 + imag_shift]
!! - the function values obtained with the pade model are compared against the 
!!   correct ones using the sum of the residuals
program test_gx_analytic_continuation

    use gx_ac, only: create_thiele_pade, evaluate_thiele_pade_at, & 
                     free_params, params

    implicit none 

    integer,      parameter :: n_parameters = 16            ! number of pade parameters (= points on imaginary axis)
    integer,      parameter :: n_points     = 100           ! number of interpolated points on real axis
    real(kind=8), parameter :: imag_shift   = 0.001d0       ! imaginary shift for interpolated points 
    logical                 :: do_greedy                    ! whether greedy algorithm is used in pade approximation
    integer                 :: precision                    ! internal floating point precision of pade fit

    type(params)                               :: params_thiele
    complex(kind=8), dimension(:), allocatable :: x_query
    complex(kind=8), dimension(:), allocatable :: y_return
    complex(kind=8), dimension(:), allocatable :: y_correct
    complex(kind=8), dimension(:), allocatable :: x_ref
    complex(kind=8), dimension(:), allocatable :: y_ref
    real(kind=8)                               :: residual_sum

    allocate(x_ref(n_parameters), y_ref(n_parameters))
    allocate(x_query(n_points), y_return(n_points), y_correct(n_points)) 

    ! get settings
    call parse_settings(do_greedy, precision)

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
                                       do_greedy=do_greedy, precision=precision)

    ! evaluate the pade interpolation model at given x points
    y_return(1:n_points) =  evaluate_thiele_pade_at(params_thiele, x_query)

    ! ----------------------------------------------------------------------------------

    ! compare the analytic continuation with the correct function values
    call evaluate_two_pole_model(n_points, x_query, y_correct)
    residual_sum = calculate_residual_sum(n_points, y_return, y_correct)
    call print_residual_sum_to_file("residual_sum.txt", residual_sum)



    ! Clean-up
    call free_params(params_thiele)
    deallocate(x_query)
    deallocate(y_return)
    deallocate(y_correct)
    deallocate(x_ref)
    deallocate(y_ref)

    contains 
        ! definition of some helper functions

        !> @brief parse the settings from the command line arguments
        !!
        !! @param[out] do_greedy -- whether greedy algorithm is used in pade fit
        !! @param[out] precision -- internal floating point precision
        subroutine parse_settings(do_greedy, precision)
            logical, intent(out) :: do_greedy 
            integer, intent(out) :: precision 

            ! internal variables 
            character(10) :: tmp 

            ! greedy or no greedy?
            call get_command_argument(1, tmp)
            select case (trim(tmp))
                case ("no_greedy")
                    do_greedy = .false.
                case ("greedy")
                    do_greedy = .true.
                case default
                    print *, "option ", trim(tmp), " not known!"
                    stop 
            end select

            ! which precision?
            call get_command_argument(2, tmp)
            read(tmp, *) precision

        end subroutine parse_settings


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


        !> @brief calculate the residual sum between the pade interpolated 
        !!        function values and the correct ones
        !!
        !! @parameter[in] n_points -- number of points
        !! @parameter[in] y        -- complex function value interpolated by pade
        !! @parameter[in] y_ref    -- correct complex function value
        !! @return        res_sum  -- sum of residuals
        real(kind=8) function calculate_residual_sum(n_points, y_pade, y_ref) result(res_sum)
            integer,                              intent(in)  :: n_points 
            complex(kind=8), dimension(n_points), intent(out) :: y_pade
            complex(kind=8), dimension(n_points), intent(out) :: y_ref

            res_sum = sum(real(abs(y_ref - y_pade), kind=8))

        end function calculate_residual_sum


        !> @brief print the residual sum to a file 
        !!
        !! @parameter[in] filename -- name of the file
        !! @parameter[in] res_sum  -- residual sum
        subroutine print_residual_sum_to_file(filename, res_sum)
            character(*), intent(in) :: filename
            real(kind=8), intent(in) :: res_sum
            
            ! internal variables
            integer :: unit_number, ios
        
            unit_number = 10
            open(unit=unit_number, file=filename, status='unknown', action='write', iostat=ios)
            write(unit_number, '(E15.8)') res_sum
            write(unit_number, *) "# options: do_greedy=", do_greedy, ", precision=", precision
            close(unit_number)

        end subroutine print_residual_sum_to_file

end program test_gx_analytic_continuation