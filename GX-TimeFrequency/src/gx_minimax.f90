!> Test Program, demonstrating how to call minimax grid routine
!>
!> Example run command: ./exe -s 32 -e ../../GX-TimeFrequency/src/ref/H2O_eigenvalues.txt -h 5 -output 32_freq_points.dat

program gx_minimax

  use iso_fortran_env, only: std_out => output_unit
  use kinds,           only: dp
  use mp2_grids,       only: gx_minimax_grid
  use minimax_exp_gw,  only: exp_gw_supported_num_points

  implicit none

!Local variables-------------------------------
!scalars
  integer :: ii, ierr, num_integ_points, unt, py_unt
  real(dp) :: emin, emax, eratio, cosft_duality_error
  character(len=1024) :: msg
  character(len=1), parameter :: ch10 = char(10)
!arrays
  real(dp) :: max_errors(3)
  real(dp), allocatable :: tau_tj(:), tau_wj(:), tj(:), wj(:)
  real(dp), allocatable :: weights_cos_tf_t_to_w(:,:), weights_cos_tf_w_to_t(:,:), weights_sin_tf_t_to_w(:,:)

! **************************************************************************************************

  emin = 1.0
  emax = 1000
  eratio = emax / emin

  open(file="gx_minimax.csv", newunit=unt, action="write")

  write(msg, "(a)") &
    "num_points ierr cosft_duality_error max_err_costf_t_to_w max_err_costf_w_to_t max_err_sintf_t_to_w eratio"
  call dumps(msg, [std_out, unt])

  do ii=1, size(exp_gw_supported_num_points)
    num_integ_points = exp_gw_supported_num_points(ii)

    call gx_minimax_grid(num_integ_points, emin, emax, &
                         tau_tj, tau_wj, tj, wj, weights_cos_tf_t_to_w, &
                         weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, max_errors, cosft_duality_error, ierr)
    write(msg, "(2(i0,2x), 5(es16.8,2x))") num_integ_points, ierr, cosft_duality_error, max_errors, eratio
    call dumps(msg, [std_out, unt])
 end do
 close(unt)

 open(file="gx_minimax.py", newunit=py_unt, action="write")
 write(py_unt, "(13a)") &
'#!/usr/bin/env python', ch10, &
'import pandas as pd', ch10, &
'import matplotlib.pyplot as plt', ch10, &
'df = pd.read_csv("gx_minimax.csv", delim_whitespace=True)', ch10, &
'ynames = ["cosft_duality_error", "max_err_costf_t_to_w",  "max_err_costf_w_to_t",  "max_err_sintf_t_to_w"] # "ierr"', ch10, &
'fig, ax_list = plt.subplots(nrows=len(ynames), ncols=1, sharex=True, squeeze=True)', ch10, &
'for ii, (yname, ax) in enumerate(zip(ynames, ax_list)):', ch10, &
'    if ii == 0:', ch10, &
'       eratio = df["eratio"][0]; max_ierr = df["ierr"].max()', ch10, &
'       fig.suptitle(f"eratio {eratio}, max_ierr: {max_ierr}")', ch10, &
'    df.plot(x="num_points", y=yname, marker="o", ax=ax, logy=True)', ch10, &
'    ax.grid(True)', ch10, &
'plt.show()', ch10
 close(py_unt)

 call execute_command_line("chmod ugo+x gx_minimax.py")
 write(std_out, "(/,a)")" Issue `./gx_minimax.py` to plot the data with matplotlib + pandas"
 call execute_command_line("./gx_minimax.py")

contains

subroutine dumps(string, units)
  character(len=*),intent(in) :: string
  integer,intent(in) :: units(:)
  integer :: jj
  do jj=1,size(units)
    write(units(jj), "(a)") trim(string)
  end do
end subroutine dumps

end program gx_minimax
