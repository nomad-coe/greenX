! **************************************************************************************************
!  Copyright (C) 2020-2025 GreenX library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************

!> @brief Program to print the maximum error of the cosine or sine transformation from frequency 
!<        to time, time to frequency and duality error

program test_gx_minimax_duality_error

    use kinds,      only: dp
    use gx_minimax, only: gx_minimax_grid, tau_npoints_supported

    implicit none

    integer                                :: n_mesh_points, ierr
    real(dp)                               :: e_transition_min, & 
                                              e_transition_max
    real(dp), dimension(:), allocatable    :: tau_mesh, tau_weights, & 
                                              freq_mesh, freq_weights
    real(dp), dimension(:, :), allocatable :: cos_tau_to_freq_weights, & 
                                              cos_freq_to_tau_weights, & 
                                              sin_tau_to_freq_weights, & 
                                              sin_freq_to_tau_weights
    real(dp), dimension(4)                 :: max_errors
    real(dp)                               :: cosft_duality_error, & 
                                              sinft_duality_error

    integer                                :: i_tau
    logical                                :: do_cosine, do_sine 
    real(dp)                               :: regularization
    real(dp), dimension(:, :), allocatable :: errors

    !Allocations
    allocate(errors(size(tau_npoints_supported),3))

    ! Parse the command line aguments 
    call parse_settings(e_transition_min, e_transition_max, do_cosine, do_sine, regularization) 

    if (do_cosine) then 
       do i_tau=1, size(tau_npoints_supported)
          n_mesh_points = tau_npoints_supported(i_tau)
          call gx_minimax_grid(n_mesh_points, e_transition_min, & 
                               e_transition_max, &
                               tau_mesh, tau_weights, &
                               freq_mesh, freq_weights, &
                               cos_tau_to_freq_weights, & 
                               cos_freq_to_tau_weights, &
                               sin_tau_to_freq_weights, &
                               max_errors(1:3), cosft_duality_error, ierr, &
                               regularization=regularization &
                                )

          ! Store the errors
          errors(i_tau,1) = max_errors(1) 
          errors(i_tau,2) = max_errors(2)
          errors(i_tau,3) = cosft_duality_error
        end do
    else if (do_sine) then
       do i_tau=1, size(tau_npoints_supported)
          n_mesh_points = tau_npoints_supported(i_tau)
          call gx_minimax_grid(n_mesh_points, e_transition_min, &
                               e_transition_max, &
                               tau_mesh, tau_weights, &
                               freq_mesh, freq_weights, &
                               cos_tau_to_freq_weights, & 
                               cos_freq_to_tau_weights, &
                               sin_tau_to_freq_weights, &
                               max_errors(1:3), cosft_duality_error, ierr, &
                               regularization=regularization, &
                               sinft_tw=sin_freq_to_tau_weights, &
                               sinft_duality_error=sinft_duality_error, & 
                               max_error_sin_wt=max_errors(4) & 
                               )
          ! Store the errors
          errors(i_tau,1) = max_errors(3)
          errors(i_tau,2) = max_errors(4)
          errors(i_tau,3) = sinft_duality_error
        end do   
    end if

    ! Print output table 
    call  print_table(e_transition_min, e_transition_max, errors)
 
    contains

    !> @brief parse the settings from the command line arguments 
    !! 
    !! @param[out] e_min -- minimun transition energy
    !! @param[out] e_max -- maximum transition energy
    !! @param[out] do_cosine -- wheter cosine dualyti error is computed
    !! @param[out] do_sine -- wheter sine duality error is computed
    !! @param[out] do_regularization -- whether Tikhonov regularization is used to compute the duality error
    !! @param[out] epsilon -- Tikhonov regularization 

    subroutine parse_settings(e_min, e_max, do_cosine, do_sine, regularization)
      real(dp), intent(inout) :: e_min, e_max, regularization
      logical,  intent(out) :: do_cosine, do_sine

      ! Local variable
      integer             :: n_args
      character(len=1024) :: tmp_val

      ! Count the number of command line arguments 
      n_args = command_argument_count()

      ! Call help function if there is not enough command line arguments
      if (n_args < 5) then
          call print_help()
      end if

      ! Read minimun transition energy
      call get_command_argument(1, tmp_val)
      if (tmp_val.eq."-emin") then
         call get_command_argument(2, tmp_val)       
         read(tmp_val, *) e_min
         if (e_min.le.0.d0) then
            write(*,*) "-emin must be greater than zero, please check"
            call print_help()
         end if 
      else
         write(*,*) "Argument at posiion 1: ", trim(tmp_val)," is not valid, calling help ..."
      end if

      ! Read maximum transition_energy
      call get_command_argument(3, tmp_val)
      if (tmp_val.eq."-emax") then
         call get_command_argument(4, tmp_val)
         read(tmp_val, *) e_max
         if (e_max.le.0.d0) then
            write(*,*) "-emax must be greater than zero, please check"
            call print_help()
         end if
      else
         write(*,*) "Argument at posiion 3: ", trim(tmp_val)," is not valid, calling help ..."
      end if
 
      call get_command_argument(5, tmp_val)       
      if (tmp_val.eq."-cosine") then
         do_cosine=.true.
         do_sine = .false.
      else if (tmp_val.eq."-sine") then
         do_cosine=.false.
         do_sine = .true.
      else 
         write(*,*) "Argument at posiion 5: ", trim(tmp_val)," is not valid, calling help ..."
            call print_help () 
      end if

      ! Set default values
      regularization = 0._dp

      ! Additional line command argument to eneble the Tikhonov regularization
      if (n_args.gt.5) then
         call get_command_argument(6, tmp_val)
         if (tmp_val.eq."-regularization") then
            if (n_args.eq.7) then 
               call get_command_argument(7, tmp_val)
               read(tmp_val, *) regularization
               if (regularization.le.0.d0) then
                  write(*,*) "-regularization must be greater than zero, please check"
                  call print_help()
               end if
            else
               regularization = 0.01_dp
            end if
         else
            write(*,*) "Argument at position 6: ", trim(tmp_val)," is not valid, calling help ..."    
            call print_help()
         end if 
      end if 

    end subroutine parse_settings

    !> @brief Print basic instructions to run the application
    subroutine print_help ()

      write(*,*) "------------------------------------------------------------------"
      write(*,*) "Usage of test_gx_minimax_duality_error.exe" 
      write(*,*) "Tabulate the cosine duality error:"
      write(*,*) "  test_gx_minimax_duality_error.exe -emin <real> -emax <real> -cosine"
      write(*,*) "Tabulate the sine duality error:"
      write(*,*) "  test_gx_minimax_duality_error.exe -emin <real> -emax <real> -sine"
      write(*,*) "If you want to use the Tikhonov regularization for cosine or sine duality error:"
      write(*,*) "  test_gx_minimax_duality_error.exe -emin <real> -emax <real> -cosine/-sine -regularization <real>"
      write(*,*) "------------------------------------------------------------------"
      write(*,*) "Example:  test_gx_minimax_duality_error.exe -emin 0.04 -emax 30.0 -cosine"
      write(*,*) "Stoping the application ..."
      stop
 
    end subroutine print_help
 
    !> @brief Print table of errors in a file
    !! 
    !! @param[out] e_min -- minimun transition energy
    !! @param[out] e_max -- maximum transition energy
    !! @param[out] errors -- freq->tau(1), tau->freq(2) and duality(3) errors

    subroutine print_table(e_min, e_max, errors)
      real(dp)                 :: e_min, e_max
      real(dp), dimension(:,:) :: errors

      ! Local variables 
      integer            :: i_tau, io_stat
      integer, parameter :: unit_number=10
      real(dp)           :: e_range 

      ! Compute the energy range
      e_range = e_max/e_min

      ! Print into a file
      open (unit_number, file='minimax_duality_error.dat', action='write', iostat=io_stat) 
      write(unit_number,'(A)') "#n   e_range           freq_to_tau_error tau_to_freq_error duality_error"
         do i_tau=1, size(tau_npoints_supported)
             n_mesh_points = tau_npoints_supported(i_tau)
             write(unit_number,'((i0,2x),4(es16.8,2x))') n_mesh_points, e_range, errors(i_tau,1:3)
         end do
      close(unit_number) 

    end subroutine print_table

end program test_gx_minimax_duality_error
