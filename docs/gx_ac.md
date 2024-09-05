---
layout: default
title: Analytic Continuation Component
tagline: Analytic Continuation
description: Analytic continuation routines for GW calculations
---

# General
This component of the GreenX library (GX-AC) implements the analytic continuation using Padé approximants.

<p align="center">
  <img src="./img/Analyticcontinuation_main.svg" alt="Performance of the GX-AC component" width="700">
</p>

Analytic continuation is a common challenge in theoretical chemistry when you have a complex analytic function $f(z)$ defined in one domain but need it in another. This issue can be addressed using Padé approximants, which are fitted to the function in one domain (such as the imaginary axis) and then evaluated in a different domain (such as the real axis). These Padé approximants are rational functions of the form
$$ f(z) \approx T_{M}(z) = \frac{A_0 + A_1z + \cdots + A_pz^p + \cdots + A_{\frac{M-1}{2}}z^{\frac{M-1}{2}}}{1 + B_1z + \cdots + B_pz^p + \cdots + B_{\frac{M}{2}}z^{\frac{M}{2}}}.$$

The GX-AC component uses the Thiele's reciprocal differences algorithm (A.B. George, Essentials of Padé Approximants, Elsevier 1975) to obtain the Padé parameters in a continued fraction form that is [equivalent](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00555) to the rational functions form above
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
This benchmark tests the numerical stability of the Padé interpolant of the GX-AC component using three model functions. In each case, a grid along the imaginary axis $x \in [0i, 1i]$  was used to determine the Padé parameters, followed by the evaluation of 1,000 function values on the real axis $x \in [0 + \eta i, 1 + \eta i]$ using the created Padé model. A small imaginary shift $\eta=0.01$ was introduced to avoid singularities in the tested pole models. The 1,000 computed points were then compared to the exact function values of the model functions to assess the mean absolute error.

<p align="center">
  <img src="./img/Analyticcontinuation_functions.svg" alt="Performance of the GX-AC component" width="1000">
</p>


### Convergence with Number of Padé Parameters 
The three model functions described above were tested with three different configurations of the GX-AC component:

- `"plain-64"`:  Thiele Padé algorithm using double-precision floating-point representation.
- `"greedy-64"`: Thiele Padé with a **greedy algorithm** and double-precision floating-point representation.
- `"plain-128"`: Thiele Padé algorithm using **quadruple-precision** (128 bit) floating points (internally).

It was observed that the greedy algorithm with quadruple precision performed similarly to "plain-128". Therefore, this configuration was left out of the plot.
The results show that using precision higher than double precision (internally) reduces the mean absolute error for all model functions. Additionally, the greedy algorithm further lowers the error compared to the plain Thiele Padé approach, at least for the tested 2-pole model.

<p align="center">
  <img src="./img/Analyticcontinuation_model_functions.svg" alt="Performance of the GX-AC component" width="700">

  <em>
  Left: Comparison of the model function with the Padé interpolated function (128 parameters) along the real axis. Right: Mean absolute error between the correct model function and interpolated functions at 1,000 test points along the real axis.
  </em>
</p>

### Performance 
Creating the Padé model (calling `create_thiele_Padé()`) scales quadratically with the number of Padé parameters (see left side of the figure below). The model settings influence the runtime as well. Using a higher precision internally will result in a higher runtime. Additionally, using the greedy algorithm for parameter evaluation will also increase the runtime compared to the plain Thiele Padé algorithm.

Evaluating the Padé model (calling `evaluate_thiele_Padé_at()`) scales linear with the number of points that are evaluated (see left side of the figure below). The type of algorithm doesn't influence the runtime but using a higher precision internally will again result in a longer runtime.

<p align="center">
  <img src="./img/Analyticcontinuation_performance.svg" alt="Performance of the GX-AC component" width="800">
</p>

## Analytic Continuation in *GW*

The [*GW* method](https://doi.org/10.3389/fchem.2019.00377) stems from many body pertubation theory and is used for calculating electronic excitations in photoemission spectroscopy. Padé approximants are used in *GW* to continue analytic functions like the [self energy](https://dx.doi.org/10.1088/1367-2630/14/5/053020) $\Sigma(\omega)$ or the [coulomb interaction](https://doi.org/10.1021/acs.jctc.3c00555) $W(\omega)$ from the imaginary to the real frequency axis. 

In this test, we present GW calculations using [FHI-aims](https://fhi-aims.org/), where either the self energy or the screened interaction is interpolated using Padé approximants from the GX-AC component. The G<sub>0</sub>W<sub>0</sub>@PBE calculations use a NAO tier 1 basis set and 400 imaginary frequency points to obtain the Padé models. For comparison, we reference a G<sub>0</sub>W<sub>0</sub>@PBE calculation using the [contour deformation](https://doi.org/10.1021/acs.jctc.8b00458) (CD) approach to accurately obtain the self-energy and screened interaction on the real axis.

The plot demonstrates that, regardless of the GX-AC component settings (greedy/non-greedy algorithm and floating-point precision), the self-energy and screened interaction can be accurately described using Padé approximants.

<p align="center">
  <img src="./img/Analyticcontinuation_gw.svg" alt="Performance of the GX-AC component" width="900">
  <em> 
  Left: Analytic continuation (AC) of the self energy of the highest occupied molecular orbital (HOMO) from the imaginary to the real frequency axis  and comparison to contour deformation.  All Padé approximants use 400 parameters.   Right: Analytic continuation of the screened interaction of a core 1s state (benzene) from the imaginary to the real frequency axis.   
  </em>
</p>


## Analytic Continuation in RT-TDDFT 





# Usage

There are two API functions that are needed in order to generate and evaluate a Thiele Padé interpolation. To create the Thiele Padé parameters call `create_thiele_Pade()` with the reference function arguments and values:
```fortran 
params_thiele = create_thiele_pade(n_par, x_ref, y_ref)
```
`x_ref`, `y_ref` must be of length `n_par`. After this step the parameters are stored in a fortran type called `params_thiele`. The parameters don't need to be accessed. In order to use the Padé model to evaluate function values with arbitrary function arguments you can use the API function `evaluate_thiele_Pade_at()`:
```fortran
y_return =  evaluate_thiele_pade_at(params_thiele, x_query)
```
`y_return` and `x_query` must be arrays of the same length. If the Padé model is not needed anymore, the parameters can be conviniently deallocated by:
```fortran 
call free_params(params_thiele)
```



### Example of a Basic Padé Interpolation

```fortran
use gx_ac, only: create_thiele_pade, evaluate_thiele_pade_at, & 
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
params_thiele = create_thiele_pade(n_par, x_ref, y_ref)

! evaluate the Padé interpolation model at given x points
y_return(1:n_fit) =  evaluate_thiele_pade_at(params_thiele, x_query)

! Clean-up
call free_params(params_thiele)
```
This is an excerpt of a stand-alone example program that can be found in `greenX/GX-AnalyticContinuation/examples/`. You can use this script to test the GX-AC component functionalities using a model function.

### Advanced usage of Padé Interpolation
By calling 
```fortran
params_thiele = create_thiele_pade(n_par, x_ref, y_ref)
y_return =  evaluate_thiele_pade_at(params_thiele, x_query)
```
the following **defaults** are used:
- thiele pade with greedy algorithm
- internal multiple precision float representation
    - turned on when GMP is linked 
    - turned off if GMP is not linked 

It is possible to change the default behavior by specifying the optional parameters `do_greedy` and `precision`. To give two examples, using the plain thiele Padé algorithm (non-greedy) with 256 bit float precision (8-fold precision): 
```fortran
params_thiele = create_thiele_pade(n_par, x_ref, y_ref, do_greedy=.false., precision=256)
y_return =  evaluate_thiele_pade_at(params_thiele, x_query)
```
or using the greedy algorithm with the faster double precision fortran implementation (doesn't make use of GMP even if it is linked) :
```fortran
params_thiele = create_thiele_pade(n_par, x_ref, y_ref, do_greedy=.true., precision=64)
y_return =  evaluate_thiele_pade_at(params_thiele, x_query)
```
All possible combinations of `do_greedy` and `precision` are supported. 

**Some considerations**:
- 64 bit precision is faster than any other precision (because only fortran is used, no GMP)
- `do_greedy=.true.` is slower than `do_greedy=.false.` 
- the routine scales $\mathcal{O}(N^2)$ in memory where $N$ is the number of Padé parameters 



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
```


<button onclick="goBack()">Go Back</button>

<script>
function goBack() {
  window.history.back();
}
</script>
