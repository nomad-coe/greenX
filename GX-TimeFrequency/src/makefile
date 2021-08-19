FC = gfortran
FFLAGS = -O3 -ffree-line-length-none -mtune=native -march=native
F90FLAGS = $(FFLAGS)
FFLAGS_NO = -O0 -ffree-line-length-none -g

PROGRAM  = get_minimax_grids

# "make" builds all
all: $(PROGRAM) realclean

get_minimax_grids: kinds.o init_types.o read_input.o init_minimax.o
	$(FC) $(FFLAGS) $(F90FLAGS) -o $@ kinds.o init_types.o read_input.o init_minimax.o  

%.o: %.f90
	$(FC) $(FFLAGS) $(F90FLAGS) -o $@ -c $<

# Utility targets
.PHONY: realclean distclean

realclean:
	rm -rf *.o *.mod

distclean: realclean
	rm -f *~ $(PROGRAM) 


