program test 

    use wigner, only: threej_table_init, threej_lookup
    use gaunt, only: r_gaunt
    use spline, only: cubic_spline
    use gauss_quadrature, only: get_gauss_legendre_grid, gauss_legendre_integrator
    use legendre_polynomial, only: evaluate_legendre_polinomial_batch, evaluate_legendre_polinomial
    use log_grid, only: create_log_grid
    use codensity_radial_function, only: calculate_codensity_radial_function
    use spherical_harmonics, only: eval_spheric_harmonic
    use constants, only: pi
    use codensity_expansion, only: beta_coefficient, capital_a_coefficient, sigma_coefficient, &
                                   codensity_expansion_coefficient, expand_codensity

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
    !call test_spherical_harmonic()

    ! test radial expansion
    !call test_radial_expansion()

    ! test codensity expansion coefficients 
    !call test_expansion_coefficients()

    ! test complete expansion
    call test_total_expansion()

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

        subroutine cartesian_to_spherical(cartesian, r, theta, phi)
            implicit none
            real(kind=8), intent(in) :: cartesian(3)  ! Input Cartesian coordinates (x, y, z)
            real(kind=8), intent(out) :: r            ! Radial distance
            real(kind=8), intent(out) :: theta        ! Polar angle
            real(kind=8), intent(out) :: phi          ! Azimuthal angle
            ! Local variables
            real(kind=8) :: x, y, z
            ! Assign Cartesian components
            x = cartesian(1)
            y = cartesian(2)
            z = cartesian(3)
            ! Compute the radial distance
            r = sqrt(x**2 + y**2 + z**2)
            ! Compute the polar angle (theta)
            if (r > 0.0d0) then
                theta = acos(z / r)
            else
                theta = 0.0d0  ! Define theta as 0 if the vector has zero length
            end if
            ! Compute the azimuthal angle (phi)
            if (x == 0.0d0 .and. y == 0.0d0) then
                phi = 0.0d0  ! Define phi as 0 if x and y are both zero
            else
                !phi = atan2(y, x)
                phi = sign(1.0d0, y)*acos(x/sqrt(x**2 + y**2))
            end if
        end subroutine cartesian_to_spherical

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

        subroutine test_radial_expansion()
            ! test if the product of the radial functions can be expandet 
            ! correctly into codensity radial functions

            integer :: l_max, i
            integer :: l, m, n_grid_points, n_int_points
            real(kind=8), dimension(3) :: r1, r2, r_eval, r_projection
            real(kind=8) :: theta_eval, phi_eval, radius_eval
            real(kind=8) :: theta_middle, phi_middle, radius_middle
            real(kind=8) :: r_a, r_b, alpha 
            real(kind=8) :: phi1, phi2
            real(kind=8) :: r_norm_grid(200)
            real(kind=8), allocatable :: expansion(:)
            real(kind=8) :: sphere_harm1, sphere_harm2, harm_sum, total_sum
            real(kind=8) :: leg_p(40), cos_theta
            real(kind=8), allocatable :: func_tmp(:), grid(:)
            type(cubic_spline) :: slater1, slater2
            type(cubic_spline), dimension(:), allocatable :: spline_g

            ! center 1
            r1 = (/-1.0d0, 0.0d0, 0.0d0/)
            ! center 2
            r2 = (/1.0d0, 0.0d0, 0.0d0/)
            ! where on the line segment should the center lie? r_c = alpha*r1 + (1-alpha)*r2
            alpha = 0.7d0

            ! other settings
            l_max = 16
            n_int_points =  50
            n_grid_points = 200
            
            ! create radial part of functions
            allocate(grid(n_grid_points), func_tmp(n_grid_points))
            call create_log_grid(40d0, n_grid_points, grid)
            call slater_function_batch(grid, 3, func_tmp)
            call slater1%create(grid, func_tmp)
            call slater_function_batch(grid, 7, func_tmp)
            call slater2%create(grid, func_tmp)

            ! project 
            r_projection = alpha*r1 + (1.0d0-alpha)*r2 
            r1 = r1 - r_projection
            r2 = r2 - r_projection

            ! length from centers to new origin (center of expansion)
            r_a = sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
            r_b = sqrt(r2(1)**2 + r2(2)**2 + r2(3)**2)

            ! mirror point of origin in line segment between a and b in spherical coordinates
            call cartesian_to_spherical(r2, radius_middle, theta_middle, phi_middle)
            
            ! get the codensity radial function 
            allocate(spline_g(l_max +1))
            call calculate_codensity_radial_function(slater1, slater2, r_a, r_b, l_max, n_int_points, spline_g)
            
            ! evaluation point in spherical coordinates
            call get_grid(-4.1d0, 4.1d0, 200, r_norm_grid)
            allocate(expansion(l_max+1))
            do i = 1, 200
                ! point where the function should be evaluated 
                r_eval =(/r_norm_grid(i), 0.0d0, 0.0d0/)
                r_eval = r_eval - r_projection
                call cartesian_to_spherical(r_eval, radius_eval, theta_eval, phi_eval)
                cos_theta = dot_product(r_eval, r2) /(r_b * radius_eval)
                ! get the expansion
                total_sum = 0.0d0
                call evaluate_legendre_polinomial(l_max, cos_theta, leg_p(1:l_max+1))
                do l = 0, l_max 
                    harm_sum = 0.0d0 
                    do m = -l, l 
                        sphere_harm1 = eval_spheric_harmonic(l, m, theta_eval, phi_eval)
                        !print *, l, m, theta_eval, phi_eval, sphere_harm1
                        sphere_harm2 = eval_spheric_harmonic(l, m, theta_middle, phi_middle)
                        harm_sum = harm_sum +  sphere_harm1*sphere_harm2
                        !print *, l, m, theta_eval, phi_eval, sphere_harm1
                    end do 
                    !harm_sum = leg_p(l+1) * (2.0d0*l + 1.0d0)
                    total_sum = total_sum + harm_sum * spline_g(l+1)%evaluate(radius_eval)
                    expansion(l+1) = total_sum * 4.0d0 * pi
                end do 
                total_sum = total_sum * 4.0d0 * pi 

                ! calculate the correct result 
                call cartesian_to_spherical(r1 - r_eval, radius_eval, theta_eval, phi_eval)
                phi1 = slater1%evaluate(radius_eval)
                call cartesian_to_spherical(r2 - r_eval, radius_eval, theta_eval, phi_eval)
                phi2 = slater2%evaluate(radius_eval)

                print *, r_norm_grid(i), phi1, phi2, phi1*phi2, expansion(:)
            end do

        end subroutine test_radial_expansion



        subroutine test_expansion_coefficients()
            real(kind=8) :: test
            real(kind=8), dimension(3) :: r1, r2
            ! center 1
            r1 = (/-1.0d0, 0.0d0, 0.0d0/)
            ! center 2
            r2 = (/1.0d0, 0.0d0, 0.0d0/)

            !test = beta_coefficient(16, 18, 5)
            !test = capital_a_coefficient(4, 3, 5, 1, 1.0d0, 1.5d0, 1.5d0, 8)
            !test = sigma_coefficient(2, 2, 0, r1, r2, 1, -1, 1, 1, 7)
            test = codensity_expansion_coefficient(3, 3, 3, 1, r1, r2, 1, -1, 1, 1, 10)

            print *, test
        end subroutine test_expansion_coefficients



        subroutine test_total_expansion()
            type(cubic_spline) :: slater1, slater2
            real(kind=8), dimension(3) :: r1, r2, r_eval
            real(kind=8), allocatable :: func_tmp(:), grid(:)
            integer :: n_grid_points
            real(kind=8) :: expansion
            real(kind=8) :: phi1, phi2, r, theta, phi

            n_grid_points = 200

            ! center 1
            r1 = (/-1.0d0, 0.0d0, 0.0d0/)
            ! center 2
            r2 = (/1.0d0, 0.0d0, 0.0d0/)
            r_eval = (/0.3d0, 0.1d0, 0.0d0/)

            ! create radial part of functions
            allocate(grid(n_grid_points), func_tmp(n_grid_points))
            call create_log_grid(40d0, n_grid_points, grid)
            call slater_function_batch(grid, 3, func_tmp)
            call slater1%create(grid, func_tmp)
            call slater_function_batch(grid, 3, func_tmp)
            call slater2%create(grid, func_tmp)

            call expand_codensity(slater1, r1, 0, 0, slater2, r2, 0, 0, 5, r_eval, expansion)
            print *, expansion

            ! correct result
            call cartesian_to_spherical(r_eval-r1, r, theta, phi)
            phi1 = slater1%evaluate(r) * eval_spheric_harmonic(0, 0, theta, phi)
            call cartesian_to_spherical(r_eval-r2, r, theta, phi)
            phi2 = slater2%evaluate(r) * eval_spheric_harmonic(0, 0, theta, phi)
            print *, phi1*phi2


        end subroutine test_total_expansion


end program test 