## GreenX library - Analytic Continuation 

The analytic continuation component provides routines to interpolate functions using the thiele pade interpolant.

> [!Note]
> **Key Features**
> - basic thiele pade algorithm
> - greedy algorithm for thiele pade to enhance numerical stability
> - arbitrary precision arithmetic using the GMP library for even more numerical stability

## Usage 

### Basic usage pade interpolation

> **Default**:
> - use the greedy algorithm
> - use 64 bit float precision (double precision) when GMP is not linked
> - use 128 bit float precision (quadrupel precision) when linked against GMP

```fortran
use gx_ac, only: create_thiele_pade, evaluate_thiele_pade_at, & 
                 free_params, params


complex(dp), allocatable :: x_ref(:), y_ref(:)

type(params)          :: params_thiele
complex(dp), dimension(:), allocatable :: x_query
complex(dp), dimension(:), allocatable :: y_return
integer                                :: n_par         ! number of pade parameters

allocate(x_ref(npar), y_ref(npar))
allocate(x_query(10), y_return(10)) 

! initialize x_ref, y_ref and x_quer

! create the interpolation model and store it in struct
params_thiele = create_thiele_pade(n_par, x_ref, y_ref)

! evaluate the interpolation model at given x points
y_return(1:10) =  evaluate_thiele_pade_at(params_thiele, x_query)

! Clean-up
call free_params(params_thiele)
```

### Advanced usage pade interpolation
e.g. to use the plain thiele pade algorithm (non-greedy) with 10 fold float precision: 
```fortran
params_thiele = create_thiele_pade(n_par, x_ref, y_ref, do_greedy=.false., precision=320)
```
e.g. to use the greedy algorithm with the faster double precision fortran implementation (doesn't make use of GMP even if it is linked) :
```fortran
params_thiele = create_thiele_pade(n_par, x_ref, y_ref, do_greedy=.true., precision=64)
```

### Availability of GMP at runtime

It is possible to check whether GMP is linked against GreenX at runtime:
```fortran
use gx_ac, only: arbitrary_precision_available, create_thiele_pade

if (arbitrary_precision_available) then
    ! this will succeed
    params_thiele = create_thiele_pade(n_par, x_ref, y_ref, do_greedy=.false., precision=320)
else if (.not. arbitrary_precision_available) then 
    ! this will result in an error
    params_thiele = create_thiele_pade(n_par, x_ref, y_ref, do_greedy=.false., precision=320)
end if 
