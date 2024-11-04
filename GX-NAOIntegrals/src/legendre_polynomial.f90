!> @brief evaluation of legendre polinomials
module legendre_polynomial

    use kinds, only: dp 
    use constants, only: pi

    implicit none 

    private
    public :: evaluate_legendre_polinomial, evaluate_legendre_polinomial_batch, &
              evaluate_renorm_assoc_leg_pol

    contains 

    !> @brief evaluate all legegendre pol. up to dregree n at point x 
    !!
    !! using the recurrence relation 
    !!        P(0,x) = 1
    !!        P(1,x) = x
    !!        P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)
    !!
    !! see https://people.sc.fsu.edu/~jburkardt/f_src/legendre_polynomial/legendre_polynomial.html
    !!
    !! @param[in]  n -- degree of legendre polynomial 
    !! @param[in]  x -- point 
    !! @param[out] p -- all legendre polynomials P_0, ..., P_n
    subroutine evaluate_legendre_polinomial(n, x, p)
        integer, intent(in) :: n 
        real(kind=dp), intent(in) :: x 
        real(kind=dp), dimension(n+1), intent(out) :: p 

        ! internal variables 
        integer :: i 

        p(1) = 1.0_dp 
        p(2) = x
        do i = 2, n 
            p(i+1) = (2.0_dp*n-1.0_dp)/dble(n) * x * p(i) - (n-1.0_dp)/dble(n) * p(i-1)
        end do 

    end subroutine evaluate_legendre_polinomial



    !> @brief evaluate all legegendre pol. up to dregree n at all points x 
    !!
    !! using the recurrence relation 
    !!        P(0,x) = 1
    !!        P(1,x) = x
    !!        P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)
    !!
    !! see https://people.sc.fsu.edu/~jburkardt/f_src/legendre_polynomial/legendre_polynomial.html
    !!
    !! @param[in]  n -- degree of legendre polynomial 
    !! @param[in]  n_points -- number of points in batch
    !! @param[in]  x -- batch of points 
    !! @param[out] p -- all legendre polynomials P_0, ..., P_n at pints x
    subroutine evaluate_legendre_polinomial_batch(n, n_points, x, p)
        integer, intent(in) :: n 
        integer, intent(in) :: n_points
        real(kind=dp), dimension(n_points), intent(in) :: x 
        real(kind=dp), dimension(n_points, n+1), intent(out) :: p 

        ! internal variables 
        integer :: i 

        p(:, 1) = 1.0_dp 
        p(:, 2) = x(:)
        do i = 2, n 
            p(:, i+1) = (2.0_dp*n-1.0_dp)/dble(n) * x(:) * p(:, i) - (n-1.0_dp)/dble(n) * p(:, i-1)
        end do 

    end subroutine evaluate_legendre_polinomial_batch 

    
    !> @brief evaluate the renormalized associated legendre polinomials
    !!
    !! renormalized \tilde{P}_l^m means that it  
    !! includes the prefactor sqrt(((2l+1)*(l-m)!) / (4pi * (l+m)!))
    !! see Numerical Recipes, 3rd edition, cambridge university press, page 294
    !!
    !! @param[in] l -- angular quantum number 
    !! @param[in] m -- magnetic quantum number (0 <= m <= l)
    !! @param[in] x -- position (-1 <= x <= 1)
    real(kind=8) function evaluate_renorm_assoc_leg_pol(l, m, x) result(ren_ass_leg_pol)
        integer, intent(in) :: l
        integer, intent(in) :: m
        real(kind=dp), intent(in) :: x 

        ! internal variables 
        integer :: i, ll 
        real(kind=dp) :: fact, oldfact, pll, pmm, pmmp1, omx2 
        
        ! consistency check 
        if ((m < 0) .or. (m > l) .or. (abs(x) > 1.0_dp)) then
            print *, "ERROR: evaluate_renorm_assoc_leg_pol: invalid input"
            stop 
        end if 
        
        ! compute \tilde{P}_m^m
        pmm = 1.0_dp 
        if (m > 0) then 
            omx2 = (1.0_dp - x)*(1.0_dp + x)
            fact = 1.0_dp 
            do i = 1, m 
                pmm  = pmm * omx2*fact/(fact + 1.0_dp)
                fact = fact + 2.0_dp 
            end do 
        end if 
        pmm = sqrt((2*m+1)*pmm / (4.0_dp * pi))
        if (mod(m, 2) == 1) then 
            pmm = -pmm
        end if 
        if (l == m) then
            ren_ass_leg_pol = pmm
            return 
        else
            ! compute \tilde{P}_{m+1}^m
            pmmp1 = x * sqrt(2.0_dp*m+3.0_dp)*pmm 
            if (l == (m+1)) then 
                ren_ass_leg_pol = pmmp1 
                return 
            else 
                ! compute \tilde{P}_{l}^m, l > m+1
                oldfact = sqrt(2.0_dp * m + 3.0_dp)
                do ll = m+2, l 
                    fact = sqrt((4.0_dp*ll*ll-1.0_dp)/dble(ll*ll-m*m))
                    pll = (x*pmmp1-pmm/oldfact)*fact 
                    oldfact = fact 
                    pmm = pmmp1
                    pmmp1 = pll
                end do 
                ren_ass_leg_pol = pll 
                return
            end if
        end if 

        print *, "ERROR: evaluate_renorm_assoc_leg_pol: something went wrong"
        stop
        
    end function evaluate_renorm_assoc_leg_pol

end module 