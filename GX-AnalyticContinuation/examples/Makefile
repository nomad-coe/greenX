# ***************************************************************************************************
#  Copyright (C) 2020-2024 GreenX library
#  This file is distributed under the terms of the APACHE2 License.
#
# ***************************************************************************************************

# this Makefile compiles the GX-AC test program pade_example.f90

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O2 -g 

# link against GX-AC component 
LIBDIR = ../../install/lib
INCLUDEDIR = ../../install/include/modules
LIBS = -L$(LIBDIR) -lgx_ac -I$(INCLUDEDIR)
RPATH = -Wl,-rpath,$(LIBDIR)

# Source file
SRC = pade_example.f90 

# Executable name
EXEC = pade_example


# Compilation rule
all: $(EXEC)

$(EXEC): $(SRC)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)  $(RPATH)

# Clean rule
clean:
	rm -f $(EXEC) *.mod *.i90 