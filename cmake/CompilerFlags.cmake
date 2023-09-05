# Set compiler flags
# CMake will append CMAKE_Fortran_FLAGS with CMAKE_Fortran_FLAGS_BUILDTYPE
# CMAKE_Fortran_FLAGS_BUILDTYPE may also have predefined values, hence initialise it 

# GCC
set(GCC_BASE
     -std='f2008'            # Fortran standard set to 2008
     -fimplicit-none         # Specify that no implicit typing is allowed
     -ffree-line-length-0    # No fixed line length
   )
set(GCC_DEBUG 
     -g               # Generate symbols
     -fbacktrace      # symbolic stack traceback
     -ffpe-trap=invalid,zero,overflow   # control over floating-point exception 
     -finit-real=nan  #  All real scalars are initialised to NaN
     -fcheck=all      # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
    )
set(GCC_RELEASE -O3)  # Level 3 optimisation. Could also consider -fpack-derived

# Intel 
set(INTEL_BASE
    -stand f08          # Fortran standard set to 2008
    -implicitnone       # Give warning when undeclared variable is used
    -diag-disable=5268  # Remove warning: Extension to standard: The text exceeds right hand column allowed on the line.
    -free               # Specify free-fortmat
   )
set(INTEL_DEBUG 
    -g           # Generate symbols
    -traceback   # symbolic stack traceback
    -fp          # Disables the ebp register in optimizations and sets the ebp register to be used as the frame pointer.
    -check all    # Checks for all runtime failures.
    -check bounds # Generates code to perform runtime checks on array subscript and character substring expressions. 
    -check-uninit #  Enables runtime checking for uninitialized variables.
    -ftrapuv      #  Set unassigned scalars as a very large integer or an invalid address
    -fpe3         # control over floating-point exception (divide by zero, overflow, invalid operation, underflow, denormalized number, positive infinity, negative infinity or a NaN)
    )

set(INTEL_RELEASE 
    -O3           # Optimsation level 3
    -no-prec-div  # Heurisitics to improves precision of floating-point divides: enables optimizations that give slightly less precise results than full IEEE division.
    -fp-model fast=2 # Semantics of floating-point calculations: Enables more aggressive optimizations on floating-point data.
    -foptimize-sibling-calls # Optimise tail recursive calls.
   )

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
   set(FF_BASE ${GCC_BASE})
   set(FF_DEBUG ${GCC_DEBUG})
   set(FF_RELEASE ${GCC_RELEASE})

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   set(FF_BASE ${INTEL_BASE})
   set(FF_DEBUG ${INTEL_DEBUG})
   set(FF_RELEASE ${INTEL_RELEASE})

else ()
     message(SEND_ERROR "flags have not been defined for this compiler: \
            ${CMAKE_Fortran_COMPILER_ID}")
endif()

string(REPLACE ";" " " FF_BASE "${FF_BASE}")
string(REPLACE ";" " " FF_DEBUG "${FF_DEBUG}")
string(REPLACE ";" " " FF_RELEASE "${FF_RELEASE}")

# Standard flags to use in all cases
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FF_BASE}")
# Initialise BUILDTYPE flags so we completely define/control
# the compiler settings
set(CMAKE_Fortran_FLAGS_DEBUG "${FF_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${FF_RELEASE}")
