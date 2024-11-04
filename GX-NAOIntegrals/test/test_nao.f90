program test 

    use wigner, only: threej_table_init, threej_lookup
    use gaunt, only: r_gaunt
    use spline, only: cubic_spline
    use gauss_quadrature, only: get_gauss_legendre_grid, gauss_legendre_integrator
    use legendre_polynomial, only: evaluate_legendre_polinomial_batch
    use log_grid, only: create_log_grid
    use codensity_radial_function, only: calculate_codensity_radial_function
    use spherical_harmonics, only: eval_spheric_harmonic
    use constants, only: pi

    implicit none 

    !integer :: i
    !integer, parameter :: n = 40
    !real(kind=8) :: a 
    !type(cubic_spline) :: my_spline
    !real(kind=8), dimension(200) :: r_grid, slater 
    !real(kind=8), dimension(122) :: r_grid_out, slater_out, y_out
    !real(kind=8), dimension(n) :: gauleg_grid, gauleg_weight, gauleg_func

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
    !call test_legendre_polinomials()

    ! test log_grid 
    !call test_log_grid()

    ! test codensity radial function
    !call test_codensity_radial_function()

    ! test spherical harmonics 
    call test_spherical_harmonic()

    contains 

        real(kind=8) function slater_function(r) result(y)
            real(kind=8), intent(in) :: r 

            ! internal variables 
            integer :: n
            real(kind=8) :: zeta
            real(kind=8) :: normc
            n = 1
            zeta = 3.0d0
            normc = (2*zeta)**n * sqrt(2*zeta/fact(2*n))
            y = normc * r**(n-1) * exp(-zeta*r)

        end function slater_function

        subroutine slater_function_batch(r, degree, y)
            real(kind=8), dimension(:), intent(in) :: r 
            integer, intent(in) :: degree
            real(kind=8), dimension(:), intent(out) :: y

            ! internal variables 
            real(kind=8) :: zeta
            real(kind=8) :: normc
            zeta = 3.0d0
            normc = (2*zeta)**degree * sqrt(2*zeta/fact(2*degree))
            y(:) = normc * r(:)**(degree-1) * exp(-zeta*r(:))

        end subroutine slater_function_batch

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



        subroutine get_grid(r_start, r_end, n, grid)
            real(kind=8), intent(in) :: r_start, r_end 
            integer, intent(in) :: n
            real(kind=8), dimension(n), intent(out) :: grid

            ! internal variables 
            real(kind=8) :: step 
            integer :: i

            step = (r_end - r_start)/dble(n-1)
            do i = 1, n 
                grid(i) = r_start + (i - 1) * step 
            end do  

        end subroutine get_grid



        real(kind=8) function fact(n)
            integer, intent(in) :: n

            integer :: i

            if (n < 0) error stop 'factorial is singular for negative integers'
            fact = 1.0d0
            do i = 2, n
                fact = fact * i
            enddo
        end function fact

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

        subroutine test_log_grid()
            integer, parameter :: n = 100
            real(kind=8), dimension(n) :: r_grid, func 
            real(kind=8), dimension(n) :: r_grid_out
            type(cubic_spline) :: my_spline
            integer :: i 

            call create_log_grid(50.0d0, n, r_grid)
            do i = 1, n 
                func(i) = slater_function(r_grid(i))
            end do
            call my_spline%create(r_grid, func)
            print *, r_grid

            call create_log_grid(40.0d0, n, r_grid_out)
            do i = 5, n 
                 print *, r_grid_out(i), slater_function(r_grid_out(i)), my_spline%evaluate(r_grid_out(i)), &
                          slater_function(r_grid_out(i)) -  my_spline%evaluate(r_grid_out(i))
            end do

        end subroutine test_log_grid

        subroutine test_codensity_radial_function()
            integer, parameter :: n = 200 
            real(kind=8), dimension(n) :: r_grid, orb1, orb2 
            type(cubic_spline) :: orb1_sp, orb2_sp 
            integer :: i
            integer, parameter :: kappa = 20
            type(cubic_spline), dimension(kappa+1) :: g_sp
            real(kind=8), dimension(200) :: norm_grid, g_out

            call create_log_grid(50.0d0, n, r_grid)
            call slater_function_batch(r_grid, 2, orb1)
            call slater_function_batch(r_grid, 4, orb2)
            call orb1_sp%create(r_grid, orb1)
            call orb2_sp%create(r_grid, orb2)

            call calculate_codensity_radial_function(orb1_sp, orb2_sp, 0.3d0, 0.7d0, kappa, 50, g_sp)

            call get_grid(0.01d0, 10.0d0, 200, norm_grid)

            do i = 1, kappa +1
                call g_sp(i)%evaluate_batch(norm_grid, g_out)
                print *, g_out
            end do


        end subroutine test_codensity_radial_function


        subroutine test_spherical_harmonic() 
            integer, parameter :: l_max = 2
            integer, parameter :: n_points = 180
            integer :: m, l, i_theta, i_phi
            real(kind=8) :: theta(n_points), phi(2*n_points), spher_harm(2*n_points, n_points) 

            theta = 0.0d0
            phi = 0.0d0

            call get_grid(0.0d0, pi, n_points, theta)
            call get_grid(0.0d0, 2.0d0*pi, 2*n_points, phi)

            l = 18
            m = 7
            do i_theta = 1, n_points
                do i_phi = 1, 2*n_points
                    spher_harm(i_phi, i_theta) = eval_spheric_harmonic(l, m, theta(i_theta), phi(i_phi))
                    !print *, l, m, theta(i_theta), phi(i_phi), spher_harm
                end do 
            end do 

            do i_theta =  1, n_points
                print *, spher_harm(:, i_theta)
            end do 
        end subroutine test_spherical_harmonic


end program test 