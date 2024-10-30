program test 

    use wigner, only: threej_table_init, threej_lookup
    use gaunt, only: r_gaunt
    use spline, only: cubic_spline
    use gauss_quadrature, only: get_gauss_legendre_grid, gauss_legendre_integrator
    use legendre_polynomial, only: evaluate_legendre_polinomial_batch

    implicit none 

    integer :: i
    integer, parameter :: n = 40
    real(kind=8) :: a 
    type(cubic_spline) :: my_spline
    real(kind=8), dimension(200) :: r_grid, slater 
    real(kind=8), dimension(122) :: r_grid_out, slater_out, y_out
    real(kind=8), dimension(n) :: gauleg_grid, gauleg_weight, gauleg_func

    !call threej_table_init(2, 40)
    !a =  threej_lookup(4, 4, 4, 0, 0, 0)

    !a = r_gaunt(10, 10, 10, 0, 0, 0)

    ! test splines
    !   call get_grid_function(0.0001d0, 10.0d0, 200, r_grid, slater)
    !   call my_spline%create(r_grid, slater)
    !   call get_grid_function(0.01d0, 9.0d0, 122, r_grid_out, slater_out)
    !   call my_spline%evaluate_batch(r_grid_out, y_out)
    !   do i = 1, 122
    !       !a = my_spline%evaluate(r_grid_out(i))
    !       print *, i, r_grid_out(i), slater_out(i), y_out(i), slater_out(i) - y_out(i)
    !   end do 


    ! test gauss integration
    !   call get_grid_function(0.0001d0, 10.0d0, 200, r_grid, slater)
    !   call my_spline%create(r_grid, slater)
    !   call get_gauss_legendre_grid(n, 0.001d0, 10.0d0, gauleg_grid, gauleg_weight)
    !   call my_spline%evaluate_batch(gauleg_grid, gauleg_func)
    !   do i = 1, n 
    !       gauleg_func(i) = (gauleg_func(i)*gauleg_grid(i))**2
    !       print *, i, gauleg_grid(i), gauleg_weight(i), gauleg_func(i)
    !   end do 
    !   a = gauss_legendre_integrator(n, gauleg_weight, gauleg_func)
    !   print *, a

    
    ! test the legendre polinomials
    call test_legendre_polinomials()


    contains 

        real(kind=8) function slater_function(r) result(y)
            real(kind=8), intent(in) :: r 

            ! internal variables 
            real(kind=8) :: zeta
            real(kind=8) :: normc
            zeta = 3.0d0
            normc = 2*zeta  * sqrt(2*zeta/2.0d0)
            y = normc * exp(-zeta*r)

        end function slater_function

        subroutine get_grid_function(r_start, r_end, n, grid, func)
            real(kind=8) :: r_start, r_end
            integer :: n 
            real(kind=8), dimension(n) :: grid, func

            real(kind=8) :: step 
            integer :: i
            step = (r_end - r_start)/(n-1)
            do i = 1, n
                grid(i) = r_start + (i-1) * step 
                func(i) = slater_function(grid(i))
            end do 

        end subroutine

        subroutine test_legendre_polinomials()
            
            integer, parameter :: n = 7 
            integer, parameter :: n_points = 300
            real(kind=8) :: p(n_points, n+1), x(n_points)

            real(kind=8), parameter :: r_start = -1.0d0
            real(kind=8), parameter :: r_end = 1.0d0

            real(kind=8) :: step 
            integer :: i, j

            step = (r_end - r_start)/(n_points-1)
            do i = 1, n_points
                x(i) = r_start + (i-1) * step 
            end do 

            call evaluate_legendre_polinomial_batch(n, n_points, x, p)

            do i = 1, n_points
                do j = 1, n+1
                    write(*, "(E20.8)", advance="no") p(i, j) 
                    write(*, "(A)", advance="no") " "
                end do 
                print *
            end do 

            
        end subroutine test_legendre_polinomials


end program test 