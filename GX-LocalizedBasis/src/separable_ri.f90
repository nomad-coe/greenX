! **************************************************************************************************
!  Copyright (C) 2020-2024 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief This module contains the subroutines for the localized basis set component of the library
! ***************************************************************************************************
module separable_ri

   use kinds,                        only: dp
   use lapack_interfaces,            only: dgemm  
   use localized_basis_types,        only: separable_ri_types
   use localized_basis_environments, only: calculate_error, power_genmat, &
                                           initialization, deallocations


   implicit none

   contains

  !> brief Compute the the three-center overlap integral (O_mn^P) using 
  !>        the separable resolution of the identity method   
  !! @param[in] n_basis_pairs: Number of orbital basis pairs (dimension 1 of ovlpXfn array)
  !! @param[in] n_loc_basbas: Number of auxiliary basis fuctions (dimension 2 of ovlp3fn array)
  !! @param[in] n_rk_points: Number of real-space grid points (dimension 2 of ovlp2fn array)
  !! param[in] ovlp_2fn: real array, the product of two NAO basis functions
  !! param[in] ovlp_3fn: real array, the three-center overlap integral over two
  !!                       NAO basis functions and one auxiliary basis function.
  !! param[out] error: real number, maximum error between Coulomb and real-space resolution of
  !!                   the identity fitting coefficients.
   subroutine gx_rirs_coefficients(n_basis_pairs,n_loc_basbas,n_rk_points, &
                                    ovlp2fn,ovlp3fn,error)

   integer n_basis_pairs,n_loc_basbas,n_rk_points

   real(kind=dp)                                        :: error
   real(kind=dp), dimension(n_basis_pairs,n_rk_points)  :: ovlp2fn
   real(kind=dp), dimension(n_basis_pairs,n_loc_basbas) :: ovlp3fn   

   ! Local variables

   type(separable_ri_types) :: ri_rs

   call initialization(ri_rs)

   ri_rs%ovlp2fn(:,:) = ovlp2fn(:,:)

   call compute_ovlp3fn(ri_rs,n_basis_pairs,n_loc_basbas,ovlp3fn)

   call calculate_error(ri_rs,n_basis_pairs,n_loc_basbas,ovlp3fn,error)

   call deallocations(ri_rs)

   end subroutine gx_rirs_coefficients

  !> brief Compute the the three-center overlap integral (O_mn^P) using 
  !>        the separable resulution of identity method  
  !! param[inout] ri_rs:  Separable resolution-of-the identity environment 
  !! param[in] ovlp_3fn:  real array, the three-center overlap integral over two
  !!                       NAO basis functions and one auxiliary basis function.
  !! @param[in] n_basis_pairs: Number of orbital basis pairs (dimension 1 of ovlp3fn array)
  !! @param[in] n_loc_basbas: Number of auxiliary basis fuctions (dimension 2 of ovlp3fn array)
   subroutine compute_ovlp3fn(ri_rs,n_basis_pairs,n_loc_basbas,ovlp3fn)
 
   type(separable_ri_types) :: ri_rs

   integer n_basis_pairs,n_loc_basbas
   real(kind=dp), dimension(n_basis_pairs, n_loc_basbas) :: ovlp3fn

   ! Compute the z coefficients 
   call get_coeff_zrs(ri_rs,n_basis_pairs,n_loc_basbas,ovlp3fn)

   ! Compute new ovlp_3fn coefficients
   call dgemm('N', 'T',n_basis_pairs,n_loc_basbas,ri_rs%n_points,1.0d0,&
              ri_rs%ovlp2fn,n_basis_pairs,ri_rs%z_coeff,n_loc_basbas,0.d0,&
              ri_rs%ovlp3fn,n_basis_pairs)

  end subroutine compute_ovlp3fn

  !> brief Compute the least-square coefficients of the separable RI method
  !>  M_Pk^k=\sum_k[(\sum_ij M_ij^P * D_ij^k') (sum_ij D_ij^k * D_ij^k')^-1]
  !! param[inout] ri_rs:  Separable resolution-of-the identity environment 
  !! param[in] ovlp_3fn:  real array, the three-center overlap integral over two
  !!                       NAO basis functions and one auxiliary basis function.
  !! @param[in] n_basis_pairs: Number of orbital basis pairs (dimension 1 of ovlp3fn array)
  !! @param[in] n_loc_basbas: Number of auxiliary basis fuctions (dimension 2 of ovlp3fn array)
   subroutine get_coeff_zrs (ri_rs,n_basis_pairs,n_loc_basbas,ovlp3fn)

   type(separable_ri_types) :: ri_rs
   
   integer n_basis_pairs,n_loc_basbas
   real(kind=dp), dimension(n_basis_pairs,n_loc_basbas) :: ovlp3fn

   !Local variables
   real(kind=dp), allocatable :: aux_mata(:,:)
   real(kind=dp), allocatable :: aux_matb(:,:)

   allocate(aux_mata(n_loc_basbas,ri_rs%n_points))
   allocate(aux_matb(ri_rs%n_points,ri_rs%n_points))

   ! Compute A = \sum_ij M_ij^P * D_ij^k'
   call dgemm('T','N',n_loc_basbas,ri_rs%n_points,n_basis_pairs,1.0d0,&
              ovlp3fn,n_basis_pairs,ri_rs%ovlp2fn,n_basis_pairs,0.0d0,&
              aux_mata,n_loc_basbas)

   ! Compute B = (sum_ij D_ij^k * D_ij^k')^-1
   call dgemm('T', 'N',ri_rs%n_points,ri_rs%n_points,n_basis_pairs,1.0d0,&
              ri_rs%ovlp2fn,n_basis_pairs,ri_rs%ovlp2fn,n_basis_pairs,0.d0,&
              aux_matb,ri_rs%n_points)

   call power_genmat(aux_matb,ri_rs%n_points,-1.d0, 1.d-10) 

   ! Compute Z = A B^-1
   call dgemm('N', 'N',n_loc_basbas,ri_rs%n_points,ri_rs%n_points,1.0d0,&
              aux_mata,n_loc_basbas,aux_matb,ri_rs%n_points,0.d0, &              
              ri_rs%z_coeff,n_loc_basbas)

   deallocate(aux_mata,aux_matb)

   end subroutine get_coeff_zrs

end module separable_ri
