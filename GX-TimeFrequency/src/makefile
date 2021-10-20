FC = gfortran
FFLAGS = -O3 -ffree-line-length-none -mtune=native -march=native
F90FLAGS = $(FFLAGS)
FFLAGS_NO = -O0 -ffree-line-length-none -g

PROGRAM = test_minimax.exe

# "make" builds all
all: $(PROGRAM) realclean

vpath %.f90 src parsers

test_minimax.exe: kinds.o parse_minimax.o parse_utils.o eigen_utilities.o minimax.o read_cmd_line.o test_minimax.o
	$(FC) $(FFLAGS) $(F90FLAGS) -o $@ kinds.o parse_minimax.o parse_utils.o eigen_utilities.o minimax.o read_cmd_line.o test_minimax.o

%.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -o $@ -c $<

# Utility targets
.PHONY: realclean distclean

realclean:
	rm -rf *.o *.mod

distclean: realclean
	rm -f *~ $(PROGRAM) 
