!> @brief evaluation of standard legendre polinomials
!!
!! using the recurrence relation 
!!        P(0,x) = 1
!!        P(1,x) = x
!!        P(n,x) = (2*n-1)/n * x * P(n-1,x) - (n-1)/n * P(n-2,x)
!!
!! see https://people.sc.fsu.edu/~jburkardt/f_src/legendre_polynomial/legendre_polynomial.html
module legendre_polynomial

    use kinds, only: dp 

    implicit none 

    private
    public :: evaluate_legendre_polinomial, evaluate_legendre_polinomial_batch

    contains 

    !> @brief evaluate all legegendre pol. up to dregree n at point x 
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

    end subroutine



    !> @brief evaluate all legegendre pol. up to dregree n at all points x 
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

    end subroutine

end module 