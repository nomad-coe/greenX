module polarizability 

  use kinds,                        only: dp
  use lapack_interfaces,            only: dgemm
  use localized_basis_environments
  use separable_ri,                 only: get_rirs_coefficients
  use gx_minimax,            only: gx_minimax_grid
  

  implicit none

  contains

  !> brief Compute the irreducible polarizability [PI(iw)]_PQ using 
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
  subroutine gx_rirs_polarizability(n_basis, n_basis_pairs, n_basbas, n_rk_points, n_states, &
                                     eigenvalues, eigenvectors, ovlp2fn, ovlp3fn, wave, error)

   integer                                              :: n_basis, n_basis_pairs, & 
                                                           n_basbas, n_states, n_rk_points 

   real(kind=dp)                                        :: error
   real(kind=dp), dimension(n_states)                   :: eigenvalues
   real(kind=dp), dimension(n_basis, n_states)          :: eigenvectors
   real(kind=dp), dimension(n_basis_pairs, n_rk_points) :: ovlp2fn
   real(kind=dp), dimension(n_basis_pairs, n_basbas)    :: ovlp3fn
   real(kind=dp), dimension(n_states, n_rk_points)      :: wave
                             

    type(polarizability_types) :: pi_pq

    ! Get the separable RI coefficients 
    call get_rirs_coefficients(pi_pq%ri_rs, n_basis_pairs, n_basbas, n_rk_points, &
                               ovlp2fn, ovlp3fn) 

    ! Get the Kohn-Sham wave function in the real-space grid
    call get_ks_wave(pi_pq, n_basis, n_states, n_rk_points, eigenvalues, & 
                     eigenvectors, wave)
                       
    ! Get minimax time-frequency grid and tranformation matrix from GreenX
    call get_minimax_grids(pi_pq)
                       
    ! Initialize space-time RPA working arrays
    call initialize_polarizability(pi_pq)

    ! Evaluate the polarizability [pi_PQ](iw)
    call evaluate_polarizability(pi_pq)    
    
    ! Deallocate working arrays
    call deallocate_polarizability(pi_pq, .true.)

    error=0.d0

  end subroutine gx_rirs_polarizability

  !> brief Compute the irreducible polarizability [PI(iw)]_PQ using 
  !>       the separable resolution of the identity method
  !! param[inout] pi_pq: polarizability environment  
  !! @param[in] n_basis: Number of orbital basis (dimension 1 of eigenvectors array) 
  !! @param[in] n_rk_points: Number of real-space grid points (dimension of wave arrays )
  !! @param[in] n_states: Number of Kohn-Sham states (dimension 1 of eingenvalues and wave
  !!                      arrays  and dimension 2 of the eigenvectors array) 
  !! param[in] eigenvalues: real array, the Kohn-Shan eigenvalues from the SCF calculation
  !! param[in] eigenvectors: real array, the Kohn-Sham eigenvectos form the SCF calculation
  !! param[in] wave: real array, the Kohn-Sham wave function in the real-space grid
  subroutine get_ks_wave(pi_pq, n_basis, n_states, n_rk_points, & 
                              eigenvalues, eigenvectors, wave) 

    type(polarizability_types) :: pi_pq

   integer                                              :: n_basis, n_states, n_rk_points    
          
    real(kind=dp), dimension(n_states)                   :: eigenvalues
    real(kind=dp), dimension(n_basis, n_states)          :: eigenvectors
    real(kind=dp), dimension(n_states, n_rk_points)      :: wave

    ! Initialization 
    pi_pq%ks%n_basis  = n_basis
    pi_pq%ks%n_states = n_states
    pi_pq%ri_rs%n_points = n_rk_points

    call initialize_kohn_sham(pi_pq)

    pi_pq%ks%eigenvalues(:,1)    = eigenvalues(:)
    pi_pq%ks%eigenvectors(:,:,1) = eigenvectors(:,:)
    pi_pq%ks%wave(:,:,1)         = wave(:,:)

  end subroutine get_ks_wave

!> brief get frequency-time minimax grids 
!  param[inout] pi_pq:  polarizability environment
! **********************************************************************
  subroutine get_minimax_grids(pi_pq)

    type(polarizability_types)               :: pi_pq


    !Local variables

    integer                                   :: info, n_points
    real(kind=8)                              :: e_tran_max, e_tran_min, duality_error
    real(kind=8), dimension(3)                :: max_errors
    real(kind=8), dimension(:),   allocatable :: freq_mesh, freq_weights
    real(kind=8), dimension(:),   allocatable :: tau_mesh, tau_weights
    real(kind=8), dimension(:,:), allocatable :: cos_tau_to_freq_weights
    real(kind=8), dimension(:,:), allocatable :: cos_freq_to_tau_weights
    real(kind=8), dimension(:,:), allocatable :: sinft_tau_to_freq_weights

    ! Initialization

    call initialize_minimax_grids(pi_pq)

    ! Get the transition energies

    ! call get_minimal_maximal_transition_energy(e_tran_min, e_tran_max)
      ! Hard coded for now
      n_points = pi_pq%minimax%n_points
      e_tran_max = 0.260
      e_tran_min = 31.725

    if (e_tran_min .le. 0.d0) then
       call register_exc("Detected metal system, not supported in minimax grid")
       return
    end if

    ! Get imaginary time and frequency points and transformations matrices

    call gx_minimax_grid(n_points, e_tran_min, e_tran_max, &
         tau_mesh, tau_weights, freq_mesh, freq_weights, &
         cos_tau_to_freq_weights, cos_freq_to_tau_weights, &
         sinft_tau_to_freq_weights, max_errors, duality_error, info)
    if (info /= 0) then
        call register_exc("Error in getting the minimax grid")
       return
    end if

    pi_pq%minimax%cos_tf(1: n_points,1: n_points) = cos_tau_to_freq_weights(1: n_points,1: n_points)
    pi_pq%minimax%omega(1: n_points)= freq_mesh(1: n_points)
    pi_pq%minimax%tau(1: n_points) = tau_mesh (1: n_points)
    pi_pq%minimax%weights(1: n_points) = freq_weights(1: n_points)

  end subroutine get_minimax_grids  

  ! ***************************************************************************
  !> brief Evaluate the polarizability [PI(iw)]_PQ in the auxiliary basis function
  !! param[inout] pi_pq: polarizability environment    
  ! ***************************************************************************
  subroutine evaluate_polarizability (pi_pq) 

    type(polarizability_types) :: pi_pq

    ! Local variables
    integer i_tau, i_omega

    ! Loop over the imaginary time grid points
    do i_tau=1, pi_pq%minimax%n_points, 1

       ! Evaluate the polarizability in the real space grid
       call evaluate_chi(pi_pq, i_tau)

       ! Transform the polarizability form real to auxiliary space
       call transform_chi_to_pi(pi_pq)

     ! Cosine transform from time to frequency domain
       do i_omega = 1, pi_pq%minimax%n_points
          pi_pq%omega(:,:,i_omega) = pi_pq%omega(:,:,i_omega) + &
                                     pi_pq%tau(:,:)*pi_pq%minimax%cos_tf(i_omega, i_tau)
       end do ! i_omega

    end do ! i_tau

  end subroutine evaluate_polarizability

  ! ***************************************************************************
  !> brief Evaluate the polarizability [chi(it)]_kk in the real-space grid
  !! param[inout] pi_pq: polarizability environment    
  ! ***************************************************************************  
  subroutine evaluate_chi(pi_pq, i_tau)

    type(polarizability_types) :: pi_pq

    integer, intent(in)        :: i_tau

   ! Local variables 

    integer                    :: i_point, i_spin, j_point
    real(kind=dp)              :: tau

    ! Initilization

    tau = pi_pq%minimax%tau(i_tau)

    ! Allocations
    allocate(pi_pq%chi%matrix(pi_pq%ri_rs%n_points,pi_pq%ri_rs%n_points))    

    do i_spin = 1, pi_pq%ks%n_spin

       call get_green_forward(pi_pq, i_spin, tau)

       call get_green_backward(pi_pq, i_spin, tau) 

    end do

    ! Get the irreducible polarizability chi_0 
    do j_point = 1, pi_pq%ri_rs%n_points
       do i_point = 1, pi_pq%ri_rs%n_points
          pi_pq%chi%matrix(i_point,j_point) = -pi_pq%chi%green_forward(i_point,j_point)* &
                                               pi_pq%chi%green_backward(i_point,j_point)
       end do
    end do

    deallocate(pi_pq%chi%green_forward,pi_pq%chi%green_backward)    

  end subroutine evaluate_chi

  ! ***************************************************************************
  !> brief Evaluate the forward green function [G(it)]_kk in the real-space grid
  !! param[inout] pi_pq: polarizability environment    
  ! *************************************************************************** 
  subroutine get_green_forward(pi_pq, i_spin, tau)

    type(polarizability_types) :: pi_pq

    integer                    :: i_spin
    real(kind=dp)              :: tau


    ! Local variables

    integer                                    :: i_state
    real(kind=dp), dimension(:,:), allocatable :: wave_occ

    ! Allocations
    allocate(pi_pq%chi%green_forward(pi_pq%ri_rs%n_points,pi_pq%ri_rs%n_points))
    allocate(wave_occ(pi_pq%ri_rs%n_points,pi_pq%ks%n_occ))
    
    
    ! Scale wave function by the eigen energy
    do i_state = 1, pi_pq%ks%n_occ
       wave_occ(:, i_state) = pi_pq%ks%wave(:,i_state,i_spin)* &
               exp(-0.5*tau*(pi_pq%ks%e_fermi - pi_pq%ks%eigenvalues(i_state,i_spin)))
    end do ! i_state

    ! Get Green function (+tau) = Psi_m(r_k)*Psi_m(r_k)
    call dgemm("n","t",pi_pq%ri_rs%n_points, pi_pq%ri_rs%n_points, pi_pq%ks%n_occ,1.d0, &
         wave_occ,pi_pq%ri_rs%n_points,wave_occ,pi_pq%ri_rs%n_points,0.d0,&
         pi_pq%chi%green_forward,pi_pq%ri_rs%n_points)

    ! Deallocations
    deallocate(wave_occ)

  end subroutine get_green_forward

  ! ***************************************************************************
  !> brief Evaluate the backward green function [G(-it)]_kk in the real-space grid
  !! param[inout] pi_pq: polarizability environment    
  ! ***************************************************************************  
  subroutine get_green_backward(pi_pq, i_spin, tau)

    type(polarizability_types) :: pi_pq
    integer                    :: i_spin
    real(kind=dp)              :: tau


    !Local variables
    integer                                   :: i_state
    real(kind=8), dimension(:,:), allocatable :: wave_virt

    ! Allocations 
    allocate(pi_pq%chi%green_backward(pi_pq%ri_rs%n_points,pi_pq%ri_rs%n_points))    
    allocate(wave_virt(pi_pq%ri_rs%n_points,pi_pq%ks%n_virt))

    ! Scale the wave function by the eigen energy
    do i_state = 1 , pi_pq%ks%n_virt
       wave_virt(:, i_state) = pi_pq%ks%wave(:,pi_pq%ks%n_occ + i_state,i_spin)* &
                exp(-0.5*tau*(pi_pq%ks%eigenvalues(pi_pq%ks%n_occ + i_state,i_spin) - & 
                pi_pq%ks%e_fermi))
    end do
    
    ! Get Green function (-tau) = Psi_a(r_k)*Psi_a(r_k')
    call dgemm("n","t",pi_pq%ri_rs%n_points, pi_pq%ri_rs%n_points, pi_pq%ks%n_virt, 1.d0, &
         wave_virt,pi_pq%ri_rs%n_points,wave_virt,pi_pq%ri_rs%n_points,0.d0,&
         pi_pq%chi%green_backward,pi_pq%ri_rs%n_points)

    deallocate (wave_virt)

  end subroutine get_green_backward  

  ! ***************************************************************************
  !> brief Transform the polarizability [chi(it)]_kk to [pi(it)]_pq 
  !! param[inout] pi_pq: polarizability environment    
  ! ***************************************************************************
  subroutine transform_chi_to_pi(pi_pq)

    type(polarizability_types) :: pi_pq

    ! Local variables
    real(kind=8), dimension(:,:), allocatable :: mat_aux

    allocate(mat_aux(pi_pq%ri_rs%n_points,pi_pq%ri_rs%basis%n_basbas))

    ! right multiplication [mat_A]_kQ = [chi_0]kk' Z_Qk'
    call dgemm("n","t",pi_pq%ri_rs%n_points, pi_pq%ri_rs%basis%n_basbas, & 
               pi_pq%ri_rs%n_points,1.d0, pi_pq%chi%matrix, pi_pq%ri_rs%n_points, &
               pi_pq%ri_rs%z_coeff, pi_pq%ri_rs%basis%n_basbas,0.d0, & 
               mat_aux, pi_pq%ri_rs%n_points)

    ! left multiplication [Pi]_PQ = Z_Pk [mat_A]_kQ
    call dgemm("n","n",pi_pq%ri_rs%basis%n_basbas, pi_pq%ri_rs%basis%n_basbas, & 
               pi_pq%ri_rs%n_points,1.d0, pi_pq%ri_rs%z_coeff, pi_pq%ri_rs%basis%n_basbas, &
               mat_aux, pi_pq%ri_rs%n_points,0.d0, pi_pq%tau,pi_pq%ri_rs%basis%n_basbas)

    deallocate(mat_aux)

  end subroutine transform_chi_to_pi
  
end module polarizability
