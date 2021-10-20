module read_cmd_line
  use kinds, only: sp, dp, long_char
  implicit none
  private

  public :: input_type

  !> Command line input parser
  type input_type
     logical                  :: have_grid_size    !! Is grid size present
     logical                  :: have_eigenvalues  !! Is eigenvalue file present
     logical                  :: have_nhomo        !! Is HOMO index present?
     logical                  :: have_kpoints      !! Is the system periodic?
     logical                  :: have_gridpath     !! Has the path to th
     logical                  :: output            !! Should grid be written to file?
     character(len=long_char) :: output_file       !! Output file for minimax
     character(len=long_char) :: eigenvaluefile    !! Input eigenvalue file 
     integer :: nhomo(2)                           !! Homo index for spin up and down channels
     integer :: npoints                            !! Number of minimax grid points
  contains
     procedure :: parse =>  parse_command_line
  end type input_type

contains

  !> Parse command line options
  subroutine parse_command_line(this)
    !> Instance of input
    class(input_type), intent(inout) :: this

    character(len=long_char) :: arg
    integer(sp) :: narg
    integer :: i
    character(len=10) :: npoints_char, nhomo_char

    this%have_grid_size = .false.
    this%have_eigenvalues = .false.
    this%have_nhomo = .false.
    this%have_kpoints = .false.
    this%output = .false.

    narg=iargc()

    if(narg < 1) then
       write(*,'(A98)') "Usage: ./exe -s gridsize -e eigenvalue_file -h nhomo -output output_dir"
       stop
    endif

    do i = 1, narg
       call getarg(i,arg)
       select case(TRIM(arg))
       case("-s")
          this%have_grid_size = .true.
          call getarg(i+1, npoints_char)
          read(npoints_char, '(I6)') this%npoints
       case("-e")
          this%have_eigenvalues = .true.
          call getarg(i+1,this%eigenvaluefile)
       case("-h")
          this%have_nhomo = .true.
          call getarg(i+1, nhomo_char)
          ! TODO(Alex) Extend to spin polarised (add assignment for nhomo(2))
          read(nhomo_char, '(I6)') this%nhomo(1)
       case("-output")
          this%output = .true.
          call getarg(i+1, this%output_file)
       end select
    enddo

    if (.not. this%have_grid_size) then
      write(*,'(A90)') "Input requires gridsize"
      stop
    endif
    if (.not. this%have_eigenvalues) then
      write(*,'(A90)') "Input requires eigenvalues file"
      stop
    endif
    if (.not. this%have_nhomo) then
      write(*,'(A90)') "Input requires HOMO index, nhomo"
      stop
    endif

  end subroutine parse_command_line

end module read_cmd_line
