## no implicit rules
.SUFFIXES:

LINALG_LIBS = -L/Users/gmatteo/local/build_gcc_12.1/lib -lopenblas -lscalapack
FCFLAGS = -g -O2 -fbacktrace -ffree-line-length-none -fallow-argument-mismatch -cpp

## definitions
FC=gfortran
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c

SOURCES=kinds.f90 mathconstants.f90 minimax_exp_k15.f90 \
	minimax_exp_k53.f90 minimax_exp.f90 minimax_exp_gw.f90 minimax_rpa.f90 mp2_grids.f90 main.f90

## main and clean targets
main: $(subst .f90,.o,$(SOURCES))
	$(FC) ${LINALG_LIBS} -o $@ $+

.PHONY: clean
clean:
	-rm -f *.o *.mod *.smod main

## compilation rules
%.o %.mod %.smod: %.f90
	$(COMPILE.f08) -o $*.o $<
	@touch $@

## .o -> .mod of the modules it uses
#main.o: one.mod

lib: main
	-rm -f libminimax.a
	ar rc libminimax.a *.o
	ranlib libminimax.a
