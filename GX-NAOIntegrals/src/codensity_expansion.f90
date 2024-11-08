!> @brief angular expansion coefficients to expand a codensity 
!!
!! subroutine naming follows the paper
!! Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
module codensity_expansion 
    
    use kinds, only: dp 
    use constants, only: pi
    use spline, only: cubic_spline
    use math, only: double_factorial, cart_to_sphere
    use gaunt, only: r_gaunt, solid_harm_prefactor, alpha_function, wigner_3j
    use spherical_harmonics, only: eval_spheric_harmonic
    use codensity_radial_function, only: calculate_codensity_radial_function

    implicit none 

    private
    public :: beta_coefficient, capital_a_coefficient, sigma_coefficient, &
              codensity_expansion_coefficient, expand_codensity

    contains 


    subroutine expand_codensity(func_1, R1, l1, m1, func_2, R2, l2, m2, l_max, R_eval, expansion)
        type(cubic_spline),          intent(in) :: func_1
        real(kind=dp), dimension(3), intent(in) :: R1 
        integer,                     intent(in) :: l1
        integer,                     intent(in) :: m1
        type(cubic_spline),          intent(in) :: func_2
        real(kind=dp), dimension(3), intent(in) :: R2
        integer,                     intent(in) :: l2
        integer,                     intent(in) :: m2
        integer,                     intent(in) :: l_max
        real(kind=dp), dimension(3), intent(in) :: R_eval
        real(kind=dp),               intent(out):: expansion

        ! internal variables 
        real(kind=dp), parameter :: alpha = 0.5_dp
        integer, parameter :: n_int_points = 50 
        integer, parameter :: n_sp_points = 200
        real(kind=dp), dimension(3) :: r_project, R1_proj, R2_proj, R_eval_proj
        integer :: i_l, i_kappa, i_sigma, i_sigma_m
        real(kind=dp) :: rad1, rad2 
        real(kind=dp) :: rad_eval, theta_eval, phi_eval
        type(cubic_spline), dimension(l_max+1) :: spline_g

        ! place the origin on the line segment between R1 and R2
        r_project = alpha * R1 + (1.0_dp-alpha) * R2
        R1_proj = R1 - r_project 
        R2_proj = R2 - r_project 
        R_eval_proj = R_eval - r_project 

        ! length from centers to new origin (center of expansion)
        rad1 = sqrt(R1_proj(1)**2 + R1_proj(2)**2 + R1_proj(3)**2)
        rad2 = sqrt(R2_proj(1)**2 + R2_proj(2)**2 + R2_proj(3)**2)
        
        ! get the expansion of the radial functions
        call calculate_codensity_radial_function(func_1, func_2, rad1, rad2, l_max, n_int_points, spline_g)

        ! convert evaluation point to spherical coordinates
        call cart_to_sphere(R_eval_proj, rad_eval, theta_eval, phi_eval)

        ! do the expansion 
        expansion = 0.0_dp 
        do i_l = 0, l_max
            do i_kappa = 0, l_max 
                do i_sigma = 0, l_max 
                    do i_sigma_m = -i_sigma, i_sigma 
                        expansion = expansion + codensity_expansion_coefficient(i_l, i_kappa, i_sigma, i_sigma_m, R1_proj, R2_proj, l1, m1, l2, m2, l_max) &
                                    * rad_eval**i_l * spline_g(i_kappa+1)%evaluate(rad_eval) &
                                    * eval_spheric_harmonic(i_sigma, i_sigma_m, theta_eval, phi_eval)
                        print *, i_l, i_kappa, i_sigma, expansion
                    end do 
                end do 
            end do 
        end do

    end subroutine expand_codensity

    !> @brief get the coefficients for the codensity expansion X_{l\kappa\Sigma\Sigma_m}
    !!
    !! equation (44) in 
    !! Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
    !!
    !!
    !! @param[in] l -- a angular momentum
    !! @param[in] kappa -- a angular momentum
    !! @param[in] sigma -- a angular momentum
    !! @param[in] sigma_m -- magnetic momentum
    !! @param[in] R1 -- coordinate of function 1
    !! @param[in] R2 -- coordinate of function 2
    !! @param[in] l1 -- angular momentum of function 1
    !! @param[in] m1 -- magnetic momentum of function 1
    !! @param[in] l2 -- angular momentum of function 2
    !! @param[in] m2 -- magnetic momentum of function 2
    !! @param[in] l_max -- maximum angular momentum in expansion
    !! @result          X coefficient 
    real(kind=dp) function codensity_expansion_coefficient(l, kappa, sigma, sigma_m, R1, R2, l1, m1, l2, m2, l_max) result(x)
        integer, intent(in) :: l
        integer, intent(in) :: kappa
        integer, intent(in) :: sigma
        integer, intent(in) :: sigma_m
        real(kind=dp), dimension(3), intent(in) :: R1
        real(kind=dp), dimension(3), intent(in) :: R2
        integer, intent(in) :: l1
        integer, intent(in) :: m1
        integer, intent(in) :: l2
        integer, intent(in) :: m2
        integer, intent(in) :: l_max

        ! internal variables 
        integer :: i_L, i_M, i_nu 
        real(kind=dp) :: tmp
        real(kind=dp) :: rad_R, theta_R, phi_R

        ! eval polar coord of R = R_a - R_b; \hat{R} = \hat{R_b}
        call cart_to_sphere(R2, rad_R, theta_R, phi_R)

        x = 0.0_dp
        do i_L = 0, l_max 
            do i_M = -i_L, i_L
                tmp = 0.0_dp
                do i_nu = -kappa, kappa 
                    tmp = tmp + r_gaunt(i_L, kappa, sigma, i_M, i_nu, sigma_m) &
                                * eval_spheric_harmonic(kappa, i_nu, theta_R, phi_R)
                end do 
                x = x + tmp * sigma_coefficient(l, i_L, i_M, R1, R2, l1, m1, l2, m2, l_max)
            end do 
        end do 
        x = x * solid_harm_prefactor(kappa)**2

    end function codensity_expansion_coefficient

    !> @brief calculate \Sigma_{l_comb L M}(R)
    !!
    !! equation (41) in 
    !! Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
    !!
    !! @param[in] l_comb -- angular momnetum = lambda1 + lambda2
    !! @param[in] L -- angular momentum
    !! @param[in] M -- magnetic momentum
    !! @param[in] R1 -- coordinate of function 1
    !! @param[in] R2 -- coordinate of function 2
    !! @param[in] l1 -- angular momentum of function 1
    !! @param[in] m1 -- magnetic momentum of function 1
    !! @param[in] l2 -- angular momentum of function 2
    !! @param[in] m2 -- magnetic momentum of function 2
    !! @param[in] l_max -- maximum angular momentum in expansion
    !! @result          Sigma coefficient 
    real(kind=dp) function sigma_coefficient(l_comb, L, M, R1, R2, l1, m1, l2, m2, l_max) result(sigma)
        integer, intent(in) :: l_comb
        integer, intent(in) :: L 
        integer, intent(in) :: M 
        real(kind=dp), dimension(3), intent(in) :: R1
        real(kind=dp), dimension(3), intent(in) :: R2
        integer, intent(in) :: l1
        integer, intent(in) :: m1
        integer, intent(in) :: l2
        integer, intent(in) :: m2
        integer, intent(in) :: l_max

        ! internal variables
        integer :: lambda_1, lambda_2, mu_1, mu_2 
        real(kind=dp) :: rad_1, rad_2, theta_1, theta_2, phi_1, phi_2 

        ! get polar coordinates of the flipped vectors
        call cart_to_sphere(-R1, rad_1, theta_1, phi_1)
        call cart_to_sphere(-R2, rad_2, theta_2, phi_2)

        ! do the summation 
        sigma = 0.0_dp
        do lambda_1 = 0, l_comb 
            lambda_2 = l_comb - lambda_1 
            do mu_1 = -lambda_1, lambda_1
                do mu_2 = -lambda_2, lambda_2
                    sigma = sigma + r_gaunt(lambda_1, lambda_2, L, mu_1, mu_2, M) &
                                    * capital_a_coefficient(l1, m1, lambda_1, mu_1, rad_1, theta_1, phi_1, l_max) &
                                    * capital_a_coefficient(l2, m2, lambda_2, mu_2, rad_2, theta_2, phi_2, l_max) 
                end do 
            end do 
        end do 

    end function sigma_coefficient



    !> @brief calculate A_{lm\lambda\mu}(R)
    !!
    !! equation (39) in
    !! Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
    !!
    !! @param[in] l -- angular momentum
    !! @param[in] m -- magnetic momentum
    !! @param[in] lambda -- angular momentum
    !! @param[in] mu -- magnetic moment
    !! @param[in] r -- radius origin to center with function Y_lm
    !! @param[in] theta -- polar angle 
    !! @param[in] phi -- azimuth
    !! @param[in] l_max -- maximum angular momentum in expansion
    !! @return    coefficient A 
    real(kind=dp) function capital_a_coefficient(l, m, lambda, mu, r, theta, phi, l_max) result(A)
        integer, intent(in) :: l
        integer, intent(in) :: m
        integer, intent(in) :: lambda
        integer, intent(in) :: mu
        real(kind=dp), intent(in) :: r
        real(kind=dp), intent(in) :: theta
        real(kind=dp), intent(in) :: phi
        integer, intent(in) :: l_max

        ! internal variables 
        integer :: i_l, i_m
        real(kind=dp) :: beta, tmp

        A = 0.0_dp
        do i_l = 0, l_max
            beta = beta_coefficient(lambda, i_l, l) * r**i_l
            tmp = 0.0_dp
            do i_m = -i_l, i_l 
                tmp = tmp + r_gaunt(l, lambda, i_l, m, mu, i_m) &
                            * eval_spheric_harmonic(i_l, i_m, theta, phi)
            end do 
            A = A + beta * tmp 
        end do 

    end function capital_a_coefficient




    !> @brief returns the coeffficient \beta(L1, L2; l)
    !!
    !! equation (13) in
    !! Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
    !!
    !! @param[in] L1 -- angular quantum number 
    !! @param[in] L2 -- angular quantum number 
    !! @param[in] l  -- angular quantum number 
    real(kind=dp) function beta_coefficient(L1, L2, l) result(beta)
        integer, intent(in) :: L1
        integer, intent(in) :: L2
        integer, intent(in) :: l 

        beta = double_factorial(2*l + 1) &
                / (double_factorial(2*L1 + 1)*double_factorial(2*L2 + 1))
        
    end function beta_coefficient


end module codensity_expansion