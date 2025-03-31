! **************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief This module contains the subroutines for the localized basis set component of the library
! ***************************************************************************************************

module w_engine

  use kinds,                        only: dp
  use lapack_interfaces,            only: dgemm
  use localized_basis_types,        only: w_engine_types
  use localized_basis_environments, only: initialize_w_engine, deallocate_w_engine, &
                                          power_genmat
  use polarizability,               only: get_rirs_polarizability
  


  implicit none

  contains

  !> brief Compute the screened Coulomb matrix (w engine) in the framework of 
  !>       the separable resolution of the identity method   
  !! @param[in] n_basis: Number of orbital basis (dimension 1 of eigenvectors array) 
  !! @param[in] n_basis_pairs: Number of orbital basis pairs (dimension 1 of ovlpXfn array)
  !! @param[in] n_loc_basbas: Number of auxiliary basis fuctions (dimension 2 of ovlp3fn array)
  !! @param[in] n_rk_points: Number of real-space grid points (dimension 2 of ovlp2fn and
  !!                         wave arrays )
  !! @param[in] n_states: Number of Kohn-Sham states (dimension 1 of eingenvalues and wave
  !!                      arrays  and dimension 2 of the eigenvectors array) 
  !! param[in] eigenvalues: real array, the Kohn-Shan eigenvalues from the SCF calculation
  !! param[in] eigenvectors: real array, the Kohn-Sham eigenvectos form the SCF calculation
  !! param[in] ovlp_2fn: real array, the product of two NAO basis functions
  !! param[in] ovlp_3fn: real array, the three-center overlap integral over two
  !!                       NAO basis functions and one auxiliary basis function.
  !! param[in] wave: real array, the Kohn-Sham wave function in the real-space grid
  !! param[out] error: real number, maximum error between Coulomb and real-space resolution of
  !!                   the identity fitting coefficients.
  subroutine gx_w_engine(n_basis, n_basis_pairs, n_basbas, n_rk_points, n_states, &
                                     eigenvalues, eigenvectors, ovlp2fn, ovlp3fn, wave, error)

    integer                                              :: n_basis, n_basis_pairs, &
                                                            n_basbas, n_states, n_rk_points

    real(kind=dp)                                        :: error
    real(kind=dp), dimension(n_states)                   :: eigenvalues
    real(kind=dp), dimension(n_basis, n_states)          :: eigenvectors
    real(kind=dp), dimension(n_basis_pairs, n_rk_points) :: ovlp2fn
    real(kind=dp), dimension(n_basis_pairs, n_basbas)    :: ovlp3fn
    real(kind=dp), dimension(n_states, n_rk_points)      :: wave

    ! Local variables
    integer                                              :: i_omega

    type(w_engine_types) :: we                              

    ! Initialize space-time RPA working arrays
    call initialize_w_engine(we)

    ! Compute the polarizability [pi(iw)]_PQ
    call get_rirs_polarizability(we%pi_pq, n_basis, n_basis_pairs, n_basbas, n_rk_points, &
                                 n_states, eigenvalues, eigenvectors, ovlp2fn, ovlp3fn, wave)

    do i_omega=1, we%pi_pq%minimax%n_points

       ! Get correlation part of the screened coulomb interaction
       call get_screened_coulomb(we, i_omega) 

    end do ! i_omega
    
       ! Deallocate working arrays
    call deallocate_w_engine(we)

    error = 0.0_dp

  end  subroutine gx_w_engine

  subroutine get_screened_coulomb(we, i_omega)

    type(w_engine_types)                       :: we

    integer                                    :: i_omega

    ! Local variables
    real(kind=dp), allocatable, dimension(:,:) :: work

    ! Compute [1-pi(iw)]_PQ
    we%omega(:, :, i_omega) = we%work(:, :) - we%pi_pq%omega(:, :, i_omega)

    ! Compute 1/[1-pi(iw)]_PQ
    call power_genmat(we%omega, we%pi_pq%ri_rs%basis%n_basbas, -1.0_dp, 1.0d-10)

    ! Compute correlation part of the screened Coulomb interaction
    we%omega(:, :, i_omega) = we%omega(:, :, i_omega) - we%work(:, :)

    ! Deallocation    
    deallocate(work)
    
  end subroutine get_screened_coulomb

end module w_engine
