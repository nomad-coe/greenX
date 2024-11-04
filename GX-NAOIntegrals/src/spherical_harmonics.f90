!> @brief tools for real spherical harmonic functions Y_{lm}(theta, phi)
module spherical_harmonics 

    use kinds, only: dp 
    use legendre_polynomial, only: evaluate_renorm_assoc_leg_pol

    implicit none 

    private 
    public :: eval_spheric_harmonic, eval_spheric_harmonic_aims
    
    contains  

        !> @brief evaluate real spherical harmonic function Y_{lm}(theta, phi)
        !!
        !! Condon–Shortley phase convention is used
        !! construction by using renormalized associated legendre polynomials
        !! where the prefactor is already included  
        !!
        !! @param[in] l -- angular momentum quantum number 
        !! @param[in] m -- magnetic quantum number 
        !! @param[in] theta --  polar angle (0 <= theta <= pi)
        !! @param[in] phi --  azimuth (0 <= theta <= 2*pi)
        real(kind=dp) function eval_spheric_harmonic(l, m, theta, phi) result(sph_harm)
            integer, intent(in) :: l 
            integer, intent(in) :: m 
            real(kind=dp), intent(in) :: theta
            real(kind=dp), intent(in) :: phi 

            ! internal variables 
            real(kind=dp) :: arg

            arg = cos(theta)

            if (m .eq. 0) then 
                sph_harm = evaluate_renorm_assoc_leg_pol(l, m, arg)
            else if (m .gt. 0) then 
                sph_harm = (-1.0_dp)**m * sqrt(2.0_dp) &
                           * evaluate_renorm_assoc_leg_pol(l, m, arg) &
                           * cos(m * phi)
            else 
                sph_harm = (-1.0_dp)**m * sqrt(2.0_dp) &
                           * evaluate_renorm_assoc_leg_pol(l, abs(m), arg) &
                           * sin(abs(m) * phi)
            end if 

        end function eval_spheric_harmonic


        !> @brief evaluate real spherical harmonic function with convention of FHI-aims
        !!
        !! using the partly Condon–Shortley phase convention as in FHI-aims 
        !!
        !! @param[in] l -- angular momentum quantum number 
        !! @param[in] m -- magnetic quantum number 
        !! @param[in] theta --  polar angle (0 <= theta <= pi)
        !! @param[in] phi --  azimuth (0 <= theta <= 2*pi)
        real(kind=dp) function eval_spheric_harmonic_aims(l, m, theta, phi) result(sph_harm)
            integer, intent(in) :: l 
            integer, intent(in) :: m 
            real(kind=dp), intent(in) :: theta
            real(kind=dp), intent(in) :: phi 

            ! internal variables 
            real(kind=dp) :: arg

            arg = cos(theta)

            if (m .eq. 0) then 
                sph_harm = evaluate_renorm_assoc_leg_pol(l, m, arg)
            else if (m .gt. 0) then 
                sph_harm = sqrt(2.0_dp) &
                           * evaluate_renorm_assoc_leg_pol(l, m, arg) &
                           * cos(m * phi)
            else 
                sph_harm = (-1.0_dp)**m * sqrt(2.0_dp) &
                           * evaluate_renorm_assoc_leg_pol(l, abs(m), arg) &
                           * sin(abs(m) * phi)
            end if 

        end function eval_spheric_harmonic_aims

end module spherical_harmonics