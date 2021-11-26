!> Minimax grid type and procedures
module minimax
  use kinds, only: dp, medium_char, long_char
  use eigen_utilities, only: get_erange_non_periodic
  use parse_minimax, only: read_freq_grid

  implicit none
  private

  public :: minimax_type

  ! Preprocessor variable for path to minimax grids
  ! NOTE: This is may be limited to 132 characters
#define __MINIMAX_PATH__ MINIMAX_PATH

  !> Minimax class
  !! TODO(Alex) Extend to spin polarised
  type minimax_type
     integer                                             :: n_spin = 1   !! Number of spin channels
     integer                                             :: n_points     !! Number of grid points
     integer                                             :: n_states     !! Number of eigenvalues
     integer, dimension(2)                               :: n_homo       !! Index of HOMO for spin up and spin down
     real(kind=dp)                                       :: e_range      !! Energy range (should be named ratio?)
     real(kind=dp), dimension(:,:), allocatable          :: eigenvalues  !! Eigenvalues
     real(kind=dp), dimension(:), allocatable            :: grid_freq         !! Minimax grid
     real(kind=dp), dimension(:), allocatable            :: weights_freq      !! Grid weights
     real(kind=dp), dimension(:), allocatable            :: grid_time         !! Minimax grid
     real(kind=dp), dimension(:), allocatable            :: weights_time      !! Grid weights
     character(len=long_char) :: file_path = __MINIMAX_PATH__ !! Path to minimax grids 
   contains
     procedure :: init => initialise_minimax             !! Initialise scalar attributes
     procedure :: deallocate => deallocate_minimax       !! Deallocate arrays
     procedure :: compute_e_range => compute_e_range     !! Compute the energy range
     procedure :: read_frequency_weights => read_frequency_weights !! Read minimax frequency and weight data
     procedure :: print => print_minimax                 !! Print
     procedure :: write => write_minimax
  end type minimax_type

contains
  
  !> Initialise minimax 
 subroutine initialise_minimax(this, n_points, n_states, n_homo)
   !> Instance of minimax_type
   class(minimax_type), intent(inout):: this
   !> Number of grid points
   integer, intent(in) :: n_points
   !> Number of eigenvalues, for a single spin channel and k-point (for now)
   integer, intent(in) :: n_states
   !> Index of the HOMO for spin-up and spin-down channels
   integer, intent(in) :: n_homo(2)

   this%n_points = n_points
   this%n_states = n_states
   this%n_homo = n_homo
   allocate(this%grid_freq(this%n_points), source=0.0_dp)
   allocate(this%grid_time(this%n_points), source=0.0_dp)
   allocate(this%weights_freq(this%n_points), source=0.0_dp)
   allocate(this%weights_time(this%n_points), source=0.0_dp)
   allocate(this%eigenvalues(this%n_states, this%n_spin))

end subroutine

  !> Deallocate minimax arrays  
  subroutine deallocate_minimax(this)
    !> Instance of minimax_type
    class(minimax_type), intent(inout) :: this
    if(allocated(this%eigenvalues)) deallocate(this%eigenvalues)
    if(allocated(this%grid_freq)) deallocate(this%grid_freq)
    if(allocated(this%grid_time)) deallocate(this%grid_time)
    if(allocated(this%weights_freq)) deallocate(this%weights_freq)
    if(allocated(this%weights_time)) deallocate(this%weights_time)
  end subroutine deallocate_minimax

  !> Wrapper for calculation of energy range requirement of the minimax grid
  !> TODO(Alex) Add error handling
  subroutine compute_e_range(this, is_periodic)
    !> Instance of minimax_type
    class(minimax_type), intent(inout) :: this
    !> Are the eigenvalues for a periodic system?
    logical, intent(in) :: is_periodic

    if (is_periodic) then
      write(*, *) 'get_erange needs to be implemented for periodic systems'
      stop
    else
      this%e_range = get_erange_non_periodic(this%eigenvalues, this%n_homo)
    endif

  end subroutine compute_e_range

  !> Wrapper for parsing the minimax frequency grid points and weights
  subroutine read_frequency_weights(this, file_name)
    !> Instance of minimax_type
    class(minimax_type), intent(inout) :: this
    !> Minimax file name prepended by full path
    character(len=*), optional, intent(in) :: file_name
    !> Number of grid points, as a string
    character(len=medium_char) :: str_npoints
    !> Local file name
    character(len=long_char) :: fname

    ! Default file name
    write(str_npoints, *) this%n_points
    fname = trim(this%file_path) // '/' // trim(adjustl(str_npoints)) // "_freq_time_points.dat"

    if (present(file_name)) then
      fname = trim(adjustl(file_name))
    endif 

    call read_freq_grid(this%n_points, trim(fname), this%e_range, &
                        this%grid_freq, this%weights_freq, this%grid_time, this%weights_time)

  end subroutine read_frequency_weights

  !> Print minimax grid to stdout 
  subroutine print_minimax(this)
    !> Instance of minimax_type
    class(minimax_type), intent(inout) :: this
    integer :: i

    write(*, *) 'Minimax Grid for frequency and time'
    write(*, *) '# Frequency grid point, freq weight, time grid point, time weight:'
    do i = 1, size(this%grid_freq)
      write(*, *) this%grid_freq(i), this%weights_freq(i), this%grid_time(i), this%weights_time(i)
    enddo

  end subroutine print_minimax

  !> Write minimax grid to file
  subroutine write_minimax(this, file_name)
    !> Instance of minimax_type
    class(minimax_type), intent(inout) :: this
    !> File name 
    character(len=*), intent(in) :: file_name
    integer :: i, unit

    open(newunit=unit, file = trim(file_name))
    write(unit, *) 'Minimax Grid'
    write(unit, *) '# Grid point, weight:'
    do i = 1, size(this%grid_freq)
      write(unit, *) this%grid_freq(i), this%weights_freq(i), this%grid_time(i), this%weights_time(i)
    enddo
    close(unit)

  end subroutine write_minimax

end module minimax
