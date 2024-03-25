!{src2tex{textfont=tt}}
!!****m* ABINIT/m_build_info
!! NAME
!!  m_build_info
!!
!! FUNCTION
!!  This module contains information about this particular version of ABINIT
!!  and its build parameters (useful for debugging).
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2024 ABINIT group (Yann Pouillon, Matteo Giantomassi)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_build_info

implicit none

! Try to prevent problems with the length of the argument
! (should not exceed the line length of 132 character).
! So, always start the string to be replaced at the beginning of a line.

! Parameters set-up by Autoconf
character(len=8),parameter :: abinit_version = &
&   "10.0.0.83-58be4d"
character(len=*),parameter :: build_target   = &
&   "x86_64_darwin18.7.0_gnu9.4"

! More info on current version
character(len=*),parameter :: version_major = &
&   "10"
character(len=*),parameter :: version_minor = &
&   "0"
character(len=*),parameter :: version_micro = &
&   "0"
character(len=*),parameter :: version_build = &
&   "20240316"

! Info on compilers. Try to prevent problems with the length of the argument
! (should not exceed the line length of 132 character).
character(len=*),parameter :: cc_info    = &
&   "gnu9.4"
character(len=*),parameter :: cxx_info   = &
&   "gnu9.4"
character(len=*),parameter :: fc_info    = &
&   "gnu9.4"
character(len=*),parameter :: cc_flags   = &
&   "-O2 -g    "
character(len=*),parameter :: cxx_flags  = &
&   "-O2 -g    "
character(len=*),parameter :: fc_flags   = &
&   "-O2 -g -ffree-line-length-none   -I/opt/local/include  -I/opt/local/include  -I/opt/local/include"
character(len=*),parameter :: fc_ldflags = &
&   "   "

! Info on optimizations
character(len=*),parameter :: with_debug_flavor = &
&   "basic"
character(len=*),parameter :: with_optim_flavor = &
&   "standard"
character(len=*),parameter :: cpu_info = &
&   "unknown_unknown"

! Info on MPI
character(len=*),parameter :: with_mpi      = &
&   "/opt/local"
character(len=*),parameter :: enable_mpi_io = &
&   "yes"

! Info on openMP
character(len=*),parameter :: enable_openmp = &
&   ""

! Info on GPU
character(len=*),parameter :: with_gpu      = &
&   ""

! Info on external dependencies
character(len=*),parameter :: linalg_flavor       = &
&   "openblas"
character(len=*),parameter :: fft_flavor          = &
&   "fftw3"
character(len=*),parameter :: with_hdf5           = &
&   "yes"
character(len=*),parameter :: with_netcdf         = &
&   "yes"
character(len=*),parameter :: with_netcdf_fortran = &
&   "yes"
character(len=*),parameter :: with_libxc          = &
&   "yes"
character(len=*),parameter :: with_wannier90      = &
&   "yes"

! Info on experimental features
character(len=*),parameter :: enable_exports = &
&   ""
character(len=*),parameter :: enable_gw_dpc = &
&   "yes"

contains  !===========================================================
!!***

!!****f* ABINIT/m_build_info/dump_config
!! NAME
!!  dump_config
!!
!! FUNCTION
!!  Reports a printout of the information stored in m_build_info,
!!  useful for error messages and debugging.
!!
!! INPUTS
!!  my_unit= Fortran unit number
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

subroutine dump_config(my_unit)

implicit none

!Arguments ------------------------------------
integer,intent(in) :: my_unit

!Local variables-------------------------------

! *********************************************************************

! TODO: things that might be added through preprocessing options, e.g.
! date and time of compilation

write(my_unit,*)
write(my_unit,'(1x,a)') repeat('+',78)
write(my_unit,*)
write(my_unit,'(a)' )' === Build Information === '
write(my_unit,'(2a)')'  Version       : ',trim(abinit_version)
write(my_unit,'(2a)')'  Build target  : ',trim(build_target)
write(my_unit,'(2a)')'  Build date    : ',trim(version_build)
write(my_unit,*)
write(my_unit,'(a)' )' === Compiler Suite === '
write(my_unit,'(2a)')'  C compiler       : ',trim(cc_info)
write(my_unit,'(2a)')'  C++ compiler     : ',trim(cxx_info)
write(my_unit,'(2a)')'  Fortran compiler : ',trim(fc_info)
write(my_unit,'(2a)')'  CFLAGS           : ',trim(cc_flags)
write(my_unit,'(2a)')'  CXXFLAGS         : ',trim(cxx_flags)
write(my_unit,'(2a)')'  FCFLAGS          : ',trim(fc_flags)
write(my_unit,'(2a)')'  FC_LDFLAGS       : ',trim(fc_ldflags)
write(my_unit,*)
write(my_unit,'(a) ')' === Optimizations === '
write(my_unit,'(2a)')'  Debug level        : ',trim(with_debug_flavor)
write(my_unit,'(2a)')'  Optimization level : ',trim(with_optim_flavor)
write(my_unit,'(2a)')'  Architecture       : ',trim(cpu_info)
write(my_unit,*)
write(my_unit,'(a) ')' === Multicore === '
write(my_unit,'(2a)')'  Parallel build : ',trim(with_mpi)
write(my_unit,'(2a)')'  Parallel I/O   : ',trim(enable_mpi_io)
write(my_unit,'(2a)')'  openMP support : ',trim(enable_openmp)
write(my_unit,'(2a)')'  GPU support    : ',trim(with_gpu)
write(my_unit,*)
write(my_unit,'(a) ')' === Connectors / Fallbacks === '
write(my_unit,'(2a)')'  LINALG flavor  : ',trim(linalg_flavor)
write(my_unit,'(2a)')'  FFT flavor     : ',trim(fft_flavor)
write(my_unit,'(2a)')'  HDF5           : ',trim(with_hdf5)
write(my_unit,'(2a)')'  NetCDF         : ',trim(with_netcdf)
write(my_unit,'(2a)')'  NetCDF Fortran : ',trim(with_netcdf_fortran)
write(my_unit,'(2a)')'  LibXC          : ',trim(with_libxc)
write(my_unit,'(2a)')'  Wannier90      : ',trim(with_wannier90)
write(my_unit,*)
write(my_unit,'(a)' )' === Experimental features === '
write(my_unit,'(2a)')'  Exports             : ',trim(enable_exports)
write(my_unit,'(2a)')'  GW double-precision : ',trim(enable_gw_dpc)
write(my_unit,*)
write(my_unit,'(1x,a)') repeat('+',78)
write(my_unit,*)

end subroutine dump_config
!!***

end module m_build_info
!!***
