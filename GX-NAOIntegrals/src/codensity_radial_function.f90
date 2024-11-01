!> @brief expansion of two radial functions
!!
!! codensity radial function G_\kappa(r) as described in
!! Talman, Int. J. of Quant. Chem., 111, 2221–2227 (2011)
module codensity_radial_function

    use kinds,            only: dp 
    use spline,           only: cubic_spline
    use legendre_polynomial, only: evaluate_legendre_polinomial_batch
    use gauss_quadrature, only: get_gauss_legendre_grid
    use log_grid, only: create_log_grid

    implicit none

    private 
    public :: calculate_codensity_radial_function

    contains

    !> @brief get all codesity radial function G_0(r), ..., G_\kappa(r)
    !!
    !! @param[in] spline_a -- spline of the first function
    !! @param[in] spline_b -- spline of the second function
    !! @param[in] r_a -- distance from origin (expansion center) to center of func. a
    !! @param[in] r_b -- distance from origin (expansion center) to center of func. b
    !! @param[in] kappa -- degree of codensity radial function
    !! @param[in] n_int_points -- number of gauss-legendre grid points
    !! @param[in] spline_g -- cod. rad. fun. for all kappa as splines
    subroutine calculate_codensity_radial_function(spline_a, spline_b, r_a, r_b, kappa, n_int_points, spline_g)
        type(cubic_spline), intent(in) :: spline_a 
        type(cubic_spline), intent(in) :: spline_b 
        real(kind=dp), intent(in) :: r_a 
        real(kind=dp), intent(in) :: r_b
        integer, intent(in) :: kappa 
        integer, intent(in) :: n_int_points
        type(cubic_spline), dimension(kappa+1), intent(out) :: spline_g

        ! internal variables
        integer :: i_gauss_grid, i_kappa
        integer, parameter :: n_sp_g = 100
        real(kind=dp), parameter :: r_outer_g = 40.0
        real(kind=dp), dimension(n_int_points) :: gauleg_grid, gauleg_weights
        real(kind=dp), dimension(:, :), allocatable :: legendre_p
        real(kind=dp), dimension(n_sp_g) :: r_grid_sp, x_a, x_b, func_a, func_b
        real(kind=dp), dimension(:, :), allocatable :: func_wab, cod_rad_fun

        ! get grid for gauss integration
        call get_gauss_legendre_grid(n_int_points, -1.0_dp, 1.0_dp, gauleg_grid, gauleg_weights)

        ! get legendre polynomials P_0(r), ..., P_kappa(r)
        allocate(legendre_p(n_int_points, kappa+1))
        call evaluate_legendre_polinomial_batch(kappa, n_int_points, gauleg_grid, legendre_p)

        ! grid for splining G_\kappa
        call create_log_grid(r_outer_g, n_sp_g, r_grid_sp)

        ! evaluate the two functions on the grids
        allocate(func_wab(n_sp_g, n_int_points))
        do i_gauss_grid = 1, n_int_points
            x_a(:) = sqrt(r_grid_sp(:)**2 + 2.0_dp*r_grid_sp(:)*r_a*gauleg_grid(i_gauss_grid) + r_a**2)
            x_b(:) = sqrt(r_grid_sp(:)**2 - 2.0_dp*r_grid_sp(:)*r_b*gauleg_grid(i_gauss_grid) + r_b**2)
            call spline_a%evaluate_batch(x_a, func_a)
            call spline_b%evaluate_batch(x_b, func_b)
            func_wab(:, i_gauss_grid) = gauleg_weights(i_gauss_grid) * func_a(:) * func_b(:)
        end do

        ! do the integration 
        allocate(cod_rad_fun(n_sp_g, kappa+1))
        call dgemm('N', 'N', n_sp_g, kappa+1, n_int_points, 1.0_dp, func_wab, n_sp_g, legendre_p, n_int_points, 0.0d0, cod_rad_fun, n_sp_g)
        cod_rad_fun = 0.5 * cod_rad_fun

        ! export as spline
        do i_kappa = 1, kappa+1 
            call spline_g(i_kappa)%create(r_grid_sp, cod_rad_fun(:, i_kappa))
        end do

        deallocate(cod_rad_fun)
        deallocate(legendre_p)
        deallocate(func_wab)

    end subroutine calculate_codensity_radial_function

end module codensity_radial_function