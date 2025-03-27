! ***************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library                                                          
!  This file is distributed under the terms of the APACHE2 License.                                 
!                                                                                                   
! ***************************************************************************************************
!> \brief This module contains the types for the localized basis set component of the library
! ***************************************************************************************************
module localized_basis_types
   implicit none

   ! Subtypes

   type basis_types
      integer                                    :: n_basbas
      integer                                    :: n_basis
      integer                                    :: n_loc_basbas
      integer                                    :: n_basis_pairs
   end type basis_types

   type grids_types
     integer                                     :: m_points
     integer,      dimension(:),     allocatable :: n_points
     real(kind=8), dimension(:,:,:), allocatable :: r_points
     real(kind=8), dimension(:,:),   allocatable :: w_points
   end type grids_types

   type species_types
      integer                                    :: n_atoms
      integer                                    :: n_species
      character, dimension(:), allocatable       :: species
      real(kind=8), dimension(:,:), allocatable  :: coords
   end type

   type minimax_types
      integer                                   :: n_points ! number of points
      real(kind=8), dimension(:,:), allocatable :: cos_tf   ! transformation
      real(kind=8), dimension(:),   allocatable :: omega    ! frequency points
      real(kind=8), dimension(:),   allocatable :: tau      ! time points
      real(kind=8), dimension(:),   allocatable :: weights  ! rpa weights
   end type minimax_types

   type kohn_sham_types
      real(kind=8), dimension(:,:,:), allocatable :: wave         ! KS wave function
      real(kind=8), dimension(:,:),   allocatable :: occupied     ! occ wave function
      real(kind=8), dimension(:,:),   allocatable :: virtual      ! virt wave function
      real(kind=8), dimension(:,:,:), allocatable :: eigenvector  ! KS eigenvectors
      real(kind=8), dimension(:,:),   allocatable :: eigenvalue   ! KS eigenvalues
   end type kohn_sham_types

   type polarizabily_types
      real(kind=8), dimension(:,:),   allocatable :: chi_kk         ! Polarizability (r_k)
      real(kind=8), dimension(:,:),   allocatable :: green_forward  ! Green function (tau)
      real(kind=8), dimension(:,:),   allocatable :: green_backward ! Green function (-tau)
      real(kind=8), dimension(:,:),   allocatable :: pi_pq_tau      ! Polarizability (tau)
      real(kind=8), dimension(:,:,:), allocatable :: pi_pq_omega    ! Polarizability (omega)
      real(kind=8), dimension(:,:),   allocatable :: rirs_coeff     ! RI-RS Coeff M_{mu k}  
   end type polarizabily_types

   ! Global types

   type separable_ri_types
     integer                                     :: n_points

     real(kind=8)                                :: error

     real(kind=8), dimension(:,:), allocatable   :: ovlp2fn
     real(kind=8), dimension(:,:), allocatable   :: ovlp3fn
     real(kind=8), dimension(:,:), allocatable   :: z_coeff

     type(basis_types)                           :: basis
     type(grids_types)                           :: grids
     type(species_types)                         :: species
   end type separable_ri_types


   type screened_coulomb_types

      real(kind=8)                                :: error
      
      type(minimax_types)                         :: minimax
      type(kohn_sham_types)                       :: ks
      type(polarizabily_types)                    :: pi_pq

      real(kind=8), dimension(:,:,:), allocatable :: matrix 

   end type screened_coulomb_types 

end module localized_basis_types
