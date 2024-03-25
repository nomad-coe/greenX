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

end module localized_basis_types
