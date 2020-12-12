# All credits go to https://aoterodelaroza.github.io/devnotes/modern-fortran-makefiles/

# no implicit rules
.SUFFIXES: 

# definitions
FC=gfortran
COMPILE.f08 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c

# sources for the program
SOURCES=main.F90 m_simple_pade.F90

# main and clean targets
main: $(subst .F90,.o,$(SOURCES))
	$(FC) -o $@ $+

.PHONY: clean
clean:
	-rm -f *.o *.mod *.smod main

# compilation rules
%.o %.mod %.smod: %.F90
	$(COMPILE.f08) -o $*.o $<
	@touch $@

# .o -> .mod of the modules it uses
main.o: m_simple_pade.mod

