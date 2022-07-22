!> Program to print minimax meshes and the associated errors.

program gx_minimax

#include "gx_common.h"

  use iso_fortran_env, only: std_out => output_unit
  use kinds,           only: dp, zero, one
  use mp2_grids,       only: gx_minimax_grid
  use minimax_exp_gw,  only: exp_gw_supported_num_points

  implicit none

!Local variables-------------------------------
!scalars
integer :: ii, it, ierr, num_integ_points, unt, py_unt, nargs, iostat, ntau = -1
real(dp) :: emin = huge(one), emax = huge(one), eratio, cosft_duality_error
character(len=1024) :: msg, arg, arg_val, iomsg, command
character(len=1), parameter :: ch10 = char(10)
logical :: do_plot = .False.
!arrays
real(dp) :: max_errors(3)
real(dp), allocatable :: tau_tj(:), tau_wj(:), tj(:), wj(:)
real(dp), allocatable :: weights_cos_tf_t_to_w(:,:), weights_cos_tf_w_to_t(:,:), weights_sin_tf_t_to_w(:,:)

! **************************************************************************************************

 nargs = command_argument_count()

 if (nargs < 1) then
   call print_help()
   stop
 end if

 call get_command_argument(1, command)

 do ii=2,nargs
   call get_command_argument(ii, arg)
   if (arg(1:1) /= "-") cycle

   select case (arg)
   case ("-emin")
     call get_command_argument(ii + 1, arg_val)
     read(arg_val, *, iostat=iostat, iomsg=iomsg) emin
     call iocheck()

   case ("-emax")
     call get_command_argument(ii + 1, arg_val)
     read(arg_val, *, iostat=iostat, iomsg=iomsg) emax
     call iocheck()

   case ("-ntau")
     call get_command_argument(ii + 1, arg_val)
     read(arg_val, *, iostat=iostat, iomsg=iomsg) ntau
     call iocheck()

   case ("-plot", "-p")
     do_plot = .True.

   case ("-h", "--help")
     write(std_out, *)"Help me!"
     stop 0

   case default
     write(std_out, *)"Invalid argument:", trim(arg)
     stop 1
   end select
 end do

 if (emin == huge(one) .or. emax == huge(one)) then
   write(std_out,*)"emin and emax arguments are mandatory"; stop 1
 end if

 if (emin <= 1e-6) then
   write(std_out,*)"emin should be greater than zero"; stop 1
 end if

 eratio = emax / emin

 select case (command)
 case ("table")

   open(file="gx_minimax_table.csv", newunit=unt, action="write")
   write(msg, "(a)") &
     "num_points ierr cosft_duality_error max_err_costf_t_to_w max_err_costf_w_to_t max_err_sintf_t_to_w eratio"
   call dumps(msg, [std_out, unt])

   do ii=1, size(exp_gw_supported_num_points)
     num_integ_points = exp_gw_supported_num_points(ii)

     call gx_minimax_grid(num_integ_points, emin, emax, &
                          tau_tj, tau_wj, tj, wj, weights_cos_tf_t_to_w, &
                          weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, max_errors, cosft_duality_error, ierr)

     write(msg, "(2(i0,2x), *(es16.8,2x))") num_integ_points, ierr, cosft_duality_error, max_errors, eratio
     call dumps(msg, [std_out, unt])
   end do
   close(unt)

   ! Write python script to plot the data.
   open(file="gx_minimax_table.py", newunit=py_unt, action="write")
   write(py_unt, "(13(a,/))") &
'#!/usr/bin/env python',  &
'import pandas as pd',  &
'import matplotlib.pyplot as plt', &
'df = pd.read_csv("gx_minimax_table.csv", delim_whitespace=True)', &
'ynames = ["cosft_duality_error", "max_err_costf_t_to_w",  "max_err_costf_w_to_t",  "max_err_sintf_t_to_w"] # "ierr"', &
'fig, ax_list = plt.subplots(nrows=len(ynames), ncols=1, sharex=True, squeeze=True)', &
'for ii, (yname, ax) in enumerate(zip(ynames, ax_list)):', &
'    if ii == 0:', &
'       eratio = df["eratio"][0]; max_ierr = df["ierr"].max()', &
'       fig.suptitle(f"eratio {eratio}, max_ierr: {max_ierr}")', &
'    df.plot(x="num_points", y=yname, marker="o", ax=ax, logy=True)', &
'    ax.grid(True)', &
'plt.show()'
   close(py_unt)

   call execute_command_line("chmod ugo+x gx_minimax_table.py")
   write(std_out, "(/,a)")" Issue `./gx_minimax_table.py` to plot the data with matplotlib + pandas"
   if (do_plot) call execute_command_line("./gx_minimax_table.py")

 case ("print")

   if (ntau == -1) then
     write(std_out, *)"Print command requires ntau"; stop 1
   end if
   do ii=1, size(exp_gw_supported_num_points)
     if (ntau == exp_gw_supported_num_points(ii)) exit
   end do
   if (ii == size(exp_gw_supported_num_points) + 1) then
     write(std_out, *)"ntau ", ntau, "is not supported,"
     write(std_out, *)"Please choose among: ", exp_gw_supported_num_points
     stop 1
   end if

   call gx_minimax_grid(ntau, emin, emax, &
                        tau_tj, tau_wj, tj, wj, weights_cos_tf_t_to_w, &
                        weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, max_errors, cosft_duality_error, ierr)

   open(file="gx_minimax_print_mesh.csv", newunit=unt, action="write")
   write(msg,"(a)") &
" mesh_index tau tau_weight omega omega_weight cosft_duality_error max_err_costf_t_to_w"// &
" max_err_costf_w_to_t max_err_sintf_t_to_w eratio"
   call dumps(msg, [std_out, unt])
   do it=1,ntau
     write(msg, "(i0, *(es12.5,2x))")it, tau_tj(it), tau_wj(it), tj(it), wj(it), cosft_duality_error, max_errors, eratio
     call dumps(msg, [std_out, unt])
   end do
   close(unt)

   open(file="gx_minimax_print_costf_t_to_w.csv", newunit=unt, action="write")
   write(msg,"(a)")" mesh_index weights_cos_tf_t_to_w"
   call dumps(msg, [std_out, unt])
   do it=1,ntau
     write(msg, "(i0, *(es12.5,2x))")it, weights_cos_tf_t_to_w(:, it)
     call dumps(msg, [std_out, unt])
   end do
   close(unt)

   ! Write python script to plot the data.
   open(file="gx_minimax_print.py", newunit=py_unt, action="write")
   write(py_unt, "(13(a,/))") &
'#!/usr/bin/env python', &
'import pandas as pd', &
'import matplotlib.pyplot as plt', &
'mesh_df = pd.read_csv("gx_minimax_print_mesh.csv", delim_whitespace=True)', &
'title = ""', &
'for k in "cosft_duality_error max_err_costf_t_to_w max_err_costf_w_to_t max_err_sintf_t_to_w eratio".split():', &
'    title += f"{k}: {mesh_df[k][0]}\n"', &
'fig, ax_list = plt.subplots(nrows=2, ncols=1, sharex=True, squeeze=True)', &
'mesh_df.plot.scatter(x="mesh_index", y="tau", s="tau_weight", ax=ax_list[0])', &
'mesh_df.plot.scatter(x="mesh_index", y="omega", s="omega_weight", ax=ax_list[1])', &
'for ax in ax_list: ax.grid(True)', &
'fig.suptitle(title, fontsize=8)', &
'plt.show()'
   close(py_unt)

   call execute_command_line("chmod ugo+x gx_minimax_print.py")
   write(std_out, "(/,a)")" Issue `./gx_minimax_print.py` to plot the data with matplotlib + pandas"
   if (do_plot) call execute_command_line("./gx_minimax_print.py")

 case default
   write(std_out, *)"Invalid command: ", trim(command)
 end select

contains

subroutine print_help()
  write(std_out, "(a,/)")"Usage: gx_minimax.x command -emin NUM -emax NUM [-ntau INT] [-plot|-p]"
  write(std_out, "(a,/)")"Usage: gx_minimax.x command -emin NUM -emax NUM [-ntau INT] [-plot|-p]"
  write(std_out, "(a,/)")"To print the mesh with 6 points and erange=emax/emin, use:"
  write(std_out, "(a,/)")"    gx_minimax.x print -ntau 6 -emin 0.4 -emax 100"
  write(std_out, "(a,/)")"To print all the meshes with erange=emax/emin, use:"
  write(std_out, "(a,/)")"    gx_minimax.x table -emin 0.4 -emax 100"
end subroutine print_help

subroutine dumps(string, units)
  character(len=*),intent(in) :: string
  integer,intent(in) :: units(:)
  integer :: jj
  do jj=1,size(units)
    write(units(jj), "(a)") trim(string)
  end do
end subroutine dumps

subroutine iocheck()
  if (iostat /= 0) then
    write(std_out, *)"Error while reading argument: ", trim(arg)
    write(std_out, *)"IOMSG: ", trim(iomsg)
    stop 1
  end if
end subroutine iocheck

end program gx_minimax
