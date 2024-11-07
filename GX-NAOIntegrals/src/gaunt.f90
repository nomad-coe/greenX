!> @brief calculation of R-Gaunt coefficients
!!
!! see Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
!! and  Homeier, H. H. H.; Steinborn, E. O. JMol Struct (THEOCHEM) 1996, 368, 31.
module gaunt

    use kinds, only: dp
    use wigner, only: threej
    use constants, only: pi

    implicit none

    private
    public :: r_gaunt, wigner_3j, solid_harm_prefactor, alpha_function

    contains

    !> @prief get the R-Gaunt coefficient
    !!
    !! @param[in] l1 -- angular quantum number 
    !! @param[in] l2 -- angular quantum number 
    !! @param[in] l3 -- angular quantum number 
    !! @param[in] m1 -- magnetic quantum number 
    !! @param[in] m2 -- magnetic quantum number 
    !! @param[in] m3 -- magnetic quantum number 
    !! @return        R-Gaunt coefficient
    real(kind=dp) function r_gaunt(l1, l2, l3, m1, m2, m3) result(coeff)
        integer, intent(in) :: l1, l2, l3, m1, m2, m3 

        ! internal variables
        integer :: sign
        integer :: l1_aux, l2_aux, l3_aux, m1_aux, m2_aux, m3_aux 
        integer :: m1_abs, m2_abs, m3_abs 
    
        l1_aux = l1
        l2_aux = l2
        l3_aux = l3
        m1_aux = m1
        m2_aux = m2
        m3_aux = m3
        m1_abs = abs(m1)
        m2_abs = abs(m2)
        m3_abs = abs(m3)

        ! create permutation to ensure m1_abs < m2_abs < m3_abs 
        call order_for_m_abs(l1_aux, l2_aux, l3_aux, m1_aux, m2_aux, m3_aux, m1_abs, m2_abs, m3_abs)

        if ((m1_abs.eq.0) .and. (m2_aux.ne.m3_aux)) then 
            coeff = 0.0_dp
        else if ((m1_abs.eq.0) .and. (m2_aux.eq.m3_aux)) then
            coeff = 4.0_dp * pi * solid_harm_prefactor(l1_aux)   &
                    * solid_harm_prefactor(l2_aux)          &
                    * solid_harm_prefactor(l3_aux)          &
                    * alpha_function(l1_aux, l2_aux, l3_aux)        &
                    * wigner_3j(l1_aux, l2_aux, l3_aux, 0, m2_aux, -m3_aux)
        else 
            if (m3_abs .ne. (m1_abs+m2_abs)) then 
                !print *, "inconsistency in R-Gaunt routine ", coeff
                coeff = 0.0_dp
                return
            end if 
            if (m1_aux+m2_aux+m3_aux .gt. 0) then 
                sign = 1
            else 
                sign = -1
            end if 
            coeff = sign * 4.0_dp * pi * 1.0_dp/sqrt(2.0_dp) &
                    * solid_harm_prefactor(l1_aux)          &
                    * solid_harm_prefactor(l2_aux)          &
                    * solid_harm_prefactor(l3_aux)          &
                    * alpha_function(l1_aux, l2_aux, l3_aux)        &
                    * wigner_3j(l1_aux, l2_aux, l3_aux, -m1_abs, -m2_abs, m3_abs)
        end if 

    end function r_gaunt
    


    !> @brief reorder quantum numbers so that m1_abs < m2_abs < m3_abs
    !!
    !! @param[inout] l -- angular quantum number
    !! @param[inout] m -- magnetic quantum number
    !! @param[inout] m_abs -- abs(m)
    subroutine order_for_m_abs(l1, l2, l3, m1, m2, m3, m1_abs, m2_abs, m3_abs)
        integer, intent(inout) :: l1, l2, l3, m1, m2, m3, m1_abs, m2_abs, m3_abs

        ! internal variables
        integer :: min_idx, mid_idx, max_idx
        integer, dimension(3):: l_tmp, m_tmp, m_abs_tmp

        if ((m1_abs .le. m2_abs) .and. (m2_abs .le. m3_abs)) return 

        l_tmp = (/l1, l2, l3/)
        m_tmp = (/m1, m2, m3/)
        m_abs_tmp = (/m1_abs, m2_abs, m3_abs/)

        max_idx = maxloc(m_abs_tmp, dim=1)
        min_idx = minloc(m_abs_tmp, dim=1)
        mid_idx = 6 - (max_idx + min_idx)

        l1 = l_tmp(min_idx)
        l2 = l_tmp(mid_idx)
        l3 = l_tmp(max_idx)

        m1 = m_tmp(min_idx)
        m2 = m_tmp(mid_idx)
        m3 = m_tmp(max_idx)

        m1_abs = m_abs_tmp(min_idx)
        m2_abs = m_abs_tmp(mid_idx)
        m3_abs = m_abs_tmp(max_idx)
    end subroutine order_for_m_abs

    !> @brief returns the prefactor of a solid harmonic function Y_lm
    !!
    !! \rho_l in the Talman 2011 paper
    !!
    !! @param[in] l -- angular quantum number
    !! @result         prefator
    real(kind=dp) function solid_harm_prefactor(l) result(prefac)
        integer, intent(in) :: l 
        prefac = sqrt((2*l + 1)/(4*pi))
    end function solid_harm_prefactor

    !> @brief get the wigner-3j symbol where all magnetic quantum n. are zero
    !!
    !! @param[in] l -- angular quantum number
    !! @return         3j symbol
    real(kind=dp) function alpha_function(l1, l2, l3) result(alpha)
        integer, intent(in) :: l1, l2, l3 
        alpha = threej(2*l1, 2*l2, 2*l3, 0, 0, 0)
    end function alpha_function

    !> @brief interface to the wigner-3j module where the input is defined as 2*l and 2*m
    !!
    !! @param[inout] l -- angular quantum number
    !! @param[inout] m -- magnetic quantum number
    !! @return         3j symbol
    real(kind=dp) function wigner_3j(l1, l2, l3, m1, m2, m3) result(symbol)
        integer :: l1, l2, l3, m1, m2, m3
        symbol = threej(2*l1, 2*l2, 2*l3, 2*m1, 2*m2, 2*m3)
    end function wigner_3j

end module gaunt