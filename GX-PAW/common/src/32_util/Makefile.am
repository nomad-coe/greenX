#
# Makefile for ABINIT                                      -*- Automake -*-
# Generated by make-makefiles-corelibs on 2024/03/16 21:07:38 +0000

#
# IMPORTANT NOTE
#
# Any manual change to this file will systematically be overwritten.
# Please modify the make-makefiles-corelibs script or its config file instead.
#

AM_CFLAGS = @ABI_CPPFLAGS@
AM_CXXFLAGS = $(sd_triqs_cxxflags)
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @sd_sys_fcflags@ @fcflags_opt_32_util@

AM_CPPFLAGS = \
	@src_02_clib_fcflags@ \
	@src_10_defs_fcflags@ \
	@src_14_hidewrite_fcflags@ \
	@src_16_hideleave_fcflags@ \
	@src_17_libtetra_ext_fcflags@ \
	@src_17_yaml_out_fcflags@ \
	@src_28_numeric_noabirule_fcflags@ \
	@src_common_fcflags@ \
	@src_core_fcflags@ \
	@src_libpaw_fcflags@ \
	@sd_fft_cppflags@ \
	@sd_fft_fcflags@ \
	@sd_linalg_cppflags@ \
	@sd_linalg_fcflags@ \
	@sd_hdf5_cppflags@ \
	@sd_netcdf_cppflags@ \
	@sd_netcdf_fcflags@ \
	@sd_netcdf_fortran_cppflags@ \
	@sd_netcdf_fortran_fcflags@ \
	@fc_mod_fcflags@

# Regular source files
lib32_util_srcs = \
	m_atomdata.F90 \
	m_dummy.F90 \
	m_exp_mat.F90 \
	m_matrix.F90 \
	m_pptools.F90 \
	m_special_funcs.F90 \
	m_symtk.F90 \
	m_pstat.F90

# Source files depending on conditionals


# Library description
noinst_LIBRARIES = lib32_util.a

lib32_util_a_SOURCES= $(lib32_util_srcs)

CLEANFILES = \
	m_atomdata_cpp.f90 \
	m_dummy_cpp.f90 \
	m_exp_mat_cpp.f90 \
	m_matrix_cpp.f90 \
	m_pptools_cpp.f90 \
	m_special_funcs_cpp.f90 \
	m_symtk_cpp.f90 \
	m_pstat_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

# Dependencies (inside the directory) of directory 32_util
#
# This file has been generated by abisrc.py.
# DO NOT edit this file. All changes will be lost.
# Use `abisrc.py makemake` to regenerate the file.

CLEANFILES += \
	m_atomdata.$(MODEXT) \
	m_dummy.$(MODEXT) \
	m_exp_mat.$(MODEXT) \
	m_matrix.$(MODEXT) \
	m_pptools.$(MODEXT) \
	m_pstat.$(MODEXT) \
	m_special_funcs.$(MODEXT) \
	m_symtk.$(MODEXT) 


EXTRA_DIST += _32_util_
