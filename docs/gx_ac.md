---
layout: default
title: Analytic Continuation Component
tagline: Analytic Continuation
description: Analytic continuation routines for GW calculations
---

# General
This component of the GreenX library implements the analytic continuation using Padé approximants.

<p align="center">
  <img src="./img/Analyticcontinuation_main.svg" alt="Performance of the GX-AC component" width="700">
</p>

Analytic continuation is a common challenge in theoretical chemistry when you have a complex analytic function $f(z)$ defined in one domain but need it in another. This issue can be addressed using Padé approximants, which are fitted to the function in one domain (such as the imaginary axis) and then evaluated in a different domain (such as the real axis). These Padé approximants are rational functions of the form
$$ f(z) \approx T_{M}(z) = \frac{A_0 + A_1z + \cdots + A_pz^p + \cdots + A_{\frac{M-1}{2}}z^{\frac{M-1}{2}}}{1 + B_1z + \cdots + B_pz^p + \cdots + B_{\frac{M}{2}}z^{\frac{M}{2}}}.$$

The GX-AC component uses the Thiele's reciprocal differences algorithm (A.B. George, Essentials of Padé Approximants, Elsevier 1975) to obtain the parameters in a continued fraction form that is [equivalent](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00555) to the rational functions form above
$$T_M(z) = \cfrac{a_1}{1+ \cfrac{a_2(z - z_1)}{\quad\ddots\quad 1+ \cfrac{a_p(z-z_{p-1})}{1+\cfrac{a_{p+1}(z-z_p)}{\quad\ddots\quad 1+a_M(z-z_{M-1})}}}}$$

where $\{z_i\}$ are a set of reference points that are used to create the Padé model.  The following relation holds for every reference point:
$$ f(z_i) = T_M(z_i)\qquad i = 1, \;\dots,\; M$$

the parameters $a_i$ are obtained by recursion:
$$a_i = g_i(z_i)\qquad i = 1, \;\dots,\; M$$
$$g_p(z_i) = \begin{dcases}  f(z_i) & p=1 \\\;\\ \frac{g_{p-1}(z_{p-1})-g_{p-1}(z_i)}{(z_i - z_{p-1})g_{p-1}(z_i)} & p>1\end{dcases}$$

Padé approximants are known to be [numerical instable](https://doi.org/10.1093/imamat/25.3.267). The GX-AC component uses two strategies to numerically stabilize the interpolation. 
First, it incorporates a [greedy algorithm for Thiele Padé approximants](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00555) that minimizes the numerical error by reordering of the reference points. Additionally it is possible to use the component with a  higher internal numerical floating point precision. This helps reducing the numerical noise caused by [catastrophic cancellation](https://doi.org/10.1145/103162.103163). Catastrophic cancellation occurs when rounding errors are amplified through the subtraction of rounded numbers, such as double-precision floating-point numbers commonly used in most programs. This is implemented using the [GNU Multiple Precision (GMP) library](https://gmplib.org/) which allows floating-point operations with customizable precision.


# Benchmarks


## Model Functions
- which model functions where used 

### Convergence with Number of Parameters 

- plot with two functions

### Performance 
Creating the Padé model (calling `create_thiele_Padé()`) scales quadratically with the number of Padé parameters (see left side of the figure below). The model settings influence the runtime as well. Using a higher precision internally will result in a higher runtime. Additionally, using the greedy algorithm for parameter evaluation will also increase the runtime copared to the plain thiele Padé algorithm.

Evaluating the Padé model (calling `evaluate_thiele_Padé_at()`) scales linear with the number of points that are evaluated (see left side of the figure below). The type of algorithm doesn't influence the runtime but using a higher precision internally will again result in a longer runtime.

<p align="center">
  <img src="./img/Analyticcontinuation_performance.svg" alt="Performance of the GX-AC component" width="700">
</p>

## Analytic Continuation GW

### Self Energy
### Coulomb Interaction

## Analytic Continuation in RT-TDDFT 





# Usage

There are two API functions that are needed in order to generate and evaluate a thiele Padé interpolation.

To create the thiele Padé parameters call `create_thiele_Padé()` with the reference function arguments and values:
```fortran 
params_thiele = create_thiele_Padé(n_par, x_ref, y_ref)
```
`x_ref`, `y_ref` must be of length `n_par`. After this step the parameters are stored in the struct called `params_thiele`.

The parameters dont need to be accessed. In order to use the Padé model to evaluate function values with arbitrary function arguments you can use the API function `evaluate_thiele_Padé_at()`:
```fortran
y_return =  evaluate_thiele_Padé_at(params_thiele, x_query)
```
If the Padé model is not needed anymore, the parameters can be conviniently deallocated by:
```fortran 
call free_params(params_thiele)
```



### Basic usage Padé interpolation

> **Defaults**:
> - use the greedy algorithm
> - use 64 bit float precision (double precision) when GMP is not linked
> - use 128 bit float precision (quadrupel precision) when linked against GMP

```fortran
use gx_ac, only: create_thiele_Padé, evaluate_thiele_Padé_at, & 
                 free_params, params

type(params)                           :: params_thiele
complex(dp), dimension(:), allocatable :: x_ref
complex(dp), dimension(:), allocatable :: y_ref
complex(dp), dimension(:), allocatable :: x_query
complex(dp), dimension(:), allocatable :: y_return
integer                                :: n_par         ! number of Padé parameters
integer                                :: n_fit         ! number of fitting points

allocate(x_ref(n_par), y_ref(n_par))
allocate(x_query(n_fit), y_return(n_fit)) 

! initialize x_ref, y_ref and x_quer
...

! create the Padé interpolation model and store it in struct
params_thiele = create_thiele_Padé(n_par, x_ref, y_ref)

! evaluate the Padé interpolation model at given x points
y_return(1:n_fit) =  evaluate_thiele_Padé_at(params_thiele, x_query)

! Clean-up
call free_params(params_thiele)
```
This is an excerpt of a stand-alone example program that can be found in `greenX/GX-AnalyticContinuation/examples/`. You can use this script to test the GX-AC component using a model function.

### Advanced usage Padé interpolation
e.g. using the plain thiele Padé algorithm (non-greedy) with 256 bit float precision (8-fold precision): 
```fortran
params_thiele = create_thiele_Padé(n_par, x_ref, y_ref, do_greedy=.false., precision=256)
```
e.g. using the greedy algorithm with the faster double precision fortran implementation (doesn't make use of GMP even if it is linked) :
```fortran
params_thiele = create_thiele_Padé(n_par, x_ref, y_ref, do_greedy=.true., precision=64)
```
All possible combinations of `do_greedy` and `precision` are supported. 

**Some considerations**:
- 64 bit precision is faster than any other precision (because only fortran is used, no GMP)
- `do_greedy=.true.` is slower than `do_greedy=.false.` 
- the routines scale $\mathcal{O}(N^2)$ in memory where $N$ is the number of Padé parameters 



### Availability of GMP at runtime

It is possible to check whether GMP is linked against GreenX at runtime:
```fortran
use gx_ac, only: arbitrary_precision_available, create_thiele_Padé

if (arbitrary_precision_available) then
    ! this will succeed
    params_thiele = create_thiele_Padé(n_par, x_ref, y_ref, do_greedy=.false., precision=320)
else if (.not. arbitrary_precision_available) then 
    ! this will result in an error
    params_thiele = create_thiele_Padé(n_par, x_ref, y_ref, do_greedy=.false., precision=320)
end if   
```


<button onclick="goBack()">Go Back</button>

<script>
function goBack() {
  window.history.back();
}
</script>
