---
layout: default
title: Analytic Continuation Component
tagline: Analytic Continuation
description: Analytic continuation Component Documentation
---

# General
This component of the GreenX library (GX-AC) implements the analytic continuation using Padé approximants.

<div style="display:flex; justify-content: center; align-items: center; padding-bottom: 20px; padding-top: 20px">
  <img src="./img/Analyticcontinuation_main.svg" width="700">
</div>

Analytic continuation (AC) is a popular mathematical technique used to extend the domain of a complex analytic (holomorphic) function $f(z)$ beyond its original region of definition. For example, in many applications, a function initially defined on the imaginary axis can be analytically continued to the real axis. Such a continuation can be performed using interpolants like Padé functions, which are fitted to the function in one domain of the complex plain (such as the imaginary axis) and then evaluated in another domain (such as the real axis). This approach is rooted in the Identity Theorem, which states that if two analytic functions match on even a small part of their domain, they must be identical on the entire domain. These Padé approximants are rational functions of the form

```math
f(z) \approx T_{M}(z) = \frac{A_0 + A_1z + \cdots + A_pz^p + \cdots + A_{\frac{M-1}{2}}z^{\frac{M-1}{2}}}{1 + B_1z + \cdots + B_pz^p + \cdots + B_{\frac{M}{2}}z^{\frac{M}{2}}}.
```

The GX-AC component uses the Thiele's reciprocal differences algorithm (A.B. George, Essentials of Padé Approximants, Elsevier 1975) to obtain the Padé parameters in a continued fraction form that is [equivalent](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00555) to the rational functions form above
```math
T_M(z) = \cfrac{a_1}{1+ \cfrac{a_2(z - z_1)}{\quad\ddots\quad 1+ \cfrac{a_p(z-z_{p-1})}{1+\cfrac{a_{p+1}(z-z_p)}{\quad\ddots\quad 1+a_M(z-z_{M-1})}}}}
```

where $\{z_i\}$ are a set of reference points that are used to create the Padé model.  The following relation holds for every reference point:
```math
 f(z_i) = T_M(z_i)\qquad i = 1, \;\dots,\; M
 ```

the parameters $a_i$ are obtained by recursion:
```math 
a_i = g_i(z_i)\qquad i = 1, \;\dots,\; M 
```
```math 
g_p(z_i) = \begin{dcases}  f(z_i) & p=1 \\\;\\ \frac{g_{p-1}(z_{p-1})-g_{p-1}(z_i)}{(z_i - z_{p-1})g_{p-1}(z_i)} & p>1\end{dcases}
```

Padé approximants are known to be [numerical unstable](https://doi.org/10.1093/imamat/25.3.267). The GX-AC component uses two strategies to numerically stabilize the interpolation. First, it incorporates a [greedy algorithm for Thiele Padé approximants](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00555) that minimizes the numerical error by reordering of the reference points. A validation of the greedy algorithm can be found in this [reference](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00555). Additionally it is possible to use the component with a  higher internal numerical floating point precision. This helps reducing the numerical noise caused by [catastrophic cancellation](https://doi.org/10.1145/103162.103163). Catastrophic cancellation occurs when rounding errors are amplified through the subtraction of rounded numbers, such as double-precision floating-point numbers commonly used in most programs. This is implemented using the [GNU Multiple Precision (GMP) library](https://gmplib.org/) which allows floating-point operations with customizable precision. Furthermore, it is possible to impose various symmetries onto the Padé model using the GX-AC component. To maximize performance, the evaluation of the Padé  model uses the [Wallis algorithm](https://numerical.recipes/book.html), which minimizes the number of divisions, an operation that is computationally expensive, especially for complex floating-point numbers and even more so for higher-precision complex numbers. The component also allows symmetry-constraint Padé models, for a full list of supported symmetries see the section [Usage](#Usage).

# Benchmarks

In this benchmark section, we first analyze the effect of various parameters by using simple model functions, providing insights into their behavior, and then demonstrate practical applications in the field of ab initio electronic structure calculations. Specifically, we showcase the performance of our library for $GW$ calculations and real-time time-dependent density functional theory (RT-TDDFT) simulations.

## Model Functions
In this section we benchmark the numerical stability of the Padé interpolant of the GX-AC component using three model functions, a 2-pole model, an 8-pole model and the cosine function. A pole in the first two function refers to a singularity of the function on the real axis of $x$, e.g. the 2-pole model has two of these singularities. In each case, a grid along the imaginary axis $z \in [0i, 1i]$  was used to determine the Padé parameters, followed by the evaluation of 1,000 function values on the real axis $z \in [0 + \eta i, 1 + \eta i]$ using the created Padé model. A small imaginary shift $\eta=0.01$ was introduced to broaden the functions, this helps avoiding arbitrarily high function values or singularities in case of the pole models. The 1,000 computed points were then compared to the exact function values of the model functions to assess the mean absolute error.

<div style="display:flex; justify-content: center; align-items: center;">
  <div>
  <img src="./img/Analyticcontinuation_functions.svg"  width="1000">
  <br>
  <div style="display: block; padding: 20px; color: gray; text-align: justify;">
    <b>Figure 1</b> The left plot depicts the location of points in the complex plain, that were used to obtain the Padé parameters (in red) and also the location of points where the analytic continuation was performed (in blue). The upper table gives the equations for the three model functions that were used in the Benchmark and the lower table includes the parameters of the 8-pole model.  
  </div>
  </div>
</div>


### Convergence with Number of Padé Parameters 
The three model functions described above were tested with four different configurations of the GX-AC component:

- `"plain-64bit"`:  Thiele Padé algorithm using double-precision floating-point representation.
- `"greedy-64bit"`: Thiele Padé with a **greedy algorithm** and double-precision floating-point representation.
- `"plain-128bit"`: Thiele Padé algorithm using **quadruple-precision** (128 bit) floating points (internally).
- `"greedy-128bit"`: Thiele Padé with a **greedy algorithm** and **quadruple-precision** floating-point representation.

We found that the performance of "greedy-128bit" is similar to "plain-128bit". Therefore, the fourth configuration is not reported in Fig. 2.

Figure 2 (left column) shows the real part of the exact model functions and their corresponding Padé approximants, calculated with 128 parameters, for the three different configurations.
The right column of Figure 2 reports the error of the AC with respect to the number of Padé  parameters.  The error is defined as the residual sum between the values obtained from the Padé model and the exact analytic reference function. 
```math
\text{MAE} = \frac{1}{N}\sum_{i=0}^{N} |f(x_i + \eta i) - T(x_i + \eta i)| 
```

Starting with the 2-pole model, the exact model is well reproduced by the Padé approximant with 128 parameters for all three AC configurations (Fig. 2(a)). The plot of the MAE (Fig. 2(b)) indicates that similar errors are achieved already with less than 10  parameters because the model is relatively simple with few features. The MAE plot also reveals that the different configurations impact the error. Compared to "plain-64bit", the "greedy-64bit" algorithm reduces the MAE by a factor 5 and the "plain-128bit" by roughly a factor of 10.

Continuing with the 8-pole model, the Padé approximants accurately reproduce all features (Fig. 2(c)). Since the model function has more complexity, we observe a stronger dependence on the number of Padé parameters compared to the 2-pole model. As shown in Fig. 2(d), the MAE decreases until reaching 50–60 parameters, after which it levels off. The "plain-128bit" setting again yields the lowest error.

Turning to the cosine function, the Padé approximant with 128 parameters visibly deviates from the model function for $\text{Re}z > 0.7$ (Fig. 2(e)). The best agreement is achieved with the "plain-128bit" setting, which is also reflected in the MAE: it is an order of magnitude smaller compared to both "plain-64bit" and "greedy-64bit".

In general, we can conclude that the AC error is primarily determined by the number of Padé parameters and can be further reduced by using more than double precision. In some cases, improvements are achieved with the greedy algorithm without the need to increase floating-point precision.




<div style="display:flex; justify-content: center; align-items: center;">
  <div style="width: 600px;">
  <img src="./img/Analyticcontinuation_model_functions.svg">
  <br>
  <div style="display: block; padding: 20px; color: gray; text-align: justify;">
  <b> Figure 2</b> Left: Comparison of the model function with the Padé interpolated function (128 parameters) along the real axis. No singularities (poles) are visible because a broadening of $\eta=0.01$ was used. Right: Mean absolute error between the correct model function and interpolated functions at 1,000 test points along the real axis.  
  </div>
  </div>
</div>

### Computational performance 
Creating the Padé model (calling `create_thiele_Padé()`) scales quadratically with the number of Padé parameters (see left side of the figure below). The model settings influence the runtime as well. Using higher internal precision increases the computational cost and leads to a significantly longer runtime, with an increase of up to two orders of magnitude. Additionally, employing the greedy algorithm for parameter evaluation further increases the runtime compared to the plain Thiele-Padé algorithm, though this increase is limited to a factor of 3-4.

Evaluating the Padé model (calling `evaluate_thiele_Padé_at()`) scales linear with the number of points that are evaluated (see left side of the figure below). The type of algorithm doesn't influence the runtime but using a higher precision internally will again increase the run time by two orders of magnitude.

<div style="display:flex; justify-content: center; align-items: center;">
  <div style="width: 800px;">
  <img src="./img/Analyticcontinuation_performance.svg" alt="Performance of the GX-AC component">
  <br>
  <div style="display: block; padding: 20px; color: gray; text-align: justify;">
    <b>Figure 3</b> Performance benchmark of the GX-AC component using the 2-pole model function. Left: Left: Runtime for obtaining the Padé parameters as a function of the number of parameters used. Right: Runtime for evaluating function values using the Padé model with 100 parameters.
  </div>
  </div>
</div>

## Analytic Continuation in *GW*

The [*GW* approach](https://doi.org/10.3389/fchem.2019.00377) in many body pertubation theory is used for calculating electronic excitations in photoemission spectroscopy. Padé approximants are used in *GW* to continue analytic functions like the [self energy](https://dx.doi.org/10.1088/1367-2630/14/5/053020) $\Sigma(\omega)$ or the [coulomb interaction](https://doi.org/10.1021/acs.jctc.3c00555) $W(\omega)$ from the imaginary to the real frequency axis. Both, $\Sigma$ and $W$ exhibit poles on the real frequencies axis. Similarly as for the three model functions, we added a small broadening parameter $i\eta$ when plotting the functions in Fig. 4.

In this test, we present $GW$ calculations using [FHI-aims](https://fhi-aims.org/) packages, which is an all-electron code based on numeric atom-centered orbitals (NAOs). The self energy or the screened interaction is interpolated using Padé approximants from the GX-AC component. The G<sub>0</sub>W<sub>0</sub> are performed on top of a preceding DFT calculations with the Perdew-Burke-Ernzerhof (PBE) function ($G_0W_0$@PBE). We used NAO basis sets of tier 1 quality and 400 imaginary frequency points to obtain the Padé models. For comparison, we reference a G<sub>0</sub>W<sub>0</sub>@PBE calculation using the [contour deformation](https://doi.org/10.1021/acs.jctc.8b00458) (CD) approach. The CD technique is more accurate than AC, as it evaluates $\Sigma$ and $W$ directly on the real frequency axis. See [The GW Compendium](https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2019.00377/full) for a comparison of different frequency integration techniques.

Figure 4 shows that, regardless of the GX-AC component settings (greedy/non-greedy algorithm and floating-point precision), the self-energy and screened interaction can be accurately described using Padé approximants. The error is dominated by the number of Padé parameters. 

<div style="display:flex; justify-content: center; align-items: center;">
  <div style="width: 900px;">
  <img src="./img/Analyticcontinuation_gw.svg">
  <br>
  <div style="display: block; padding: 20px; color: gray; text-align: justify;">
  <b>Figure 4</b> Left: Analytic continuation (AC) of the self energy of the highest occupied molecular orbital (HOMO) from the imaginary to the real frequency axis  and comparison to contour deformation.  All Padé approximants use 400 parameters.   Right: Analytic continuation of the screened interaction of a core 1s state (benzene) from the imaginary to the real frequency axis.   
  </div>
  </div>
</div>


## Analytic Continuation in RT-TDDFT 

[RT-TDDFT](https://doi.org/10.1021/acs.chemrev.0c00223) is one of the most popular and computationally efficient methods to simulate electron dynamics. RT-TDDFT relies on the propagation of the electron density in the time-domain under external perturbation. The response properties, such as the dynamic polarizability tensor $\alpha^{\textnormal{el, RT}}(t)$, can be used to simulate the absorption and/or resonance raman spectroscopies. $\alpha^{\textnormal{el, RT}}(t)$ is computed from the induced dipole moment as the difference between the dipole moment of the perturbed system at RT-TDDFT time step $t$ and that of the unperturbed system. A single Fourier transformation of $\alpha^{\textnormal{el, RT}}(t)$ yields $\alpha^{\textnormal{el, RT}}(\omega)$, from which we can compute the absorption spectrum  $S(\omega)$ according to

```math 
S(\omega) = \frac{4\pi\omega}{3c}\textnormal{Tr}\left\{ \textnormal{Im}(\alpha_{\alpha\beta}^{\textnormal{el, RT}}(\omega ))\right\}  
```
where $c$ denotes the speed of light. 

The resolution of $S(\omega)$ depends on the length of the RT-TDDFT trajectory and increases with longer simulation times, as shown in Figure 5(a). RT-TDDFT calculations are therefore computationally demanding because obtaining a converged spectrum often requires long RT-TDDFT trajectories. [The use of Padé approximants](https://doi.org/10.1063/1.5051250) enables significantly shorter simulation times and higher resolution in the frequency domain.

 We demonstrate this for the naphtalene molecule. We computed the absorption spectrum computed via RT-TDDFT using the [CP2K](https://www.cp2k.org/) program package. We employed the PBE functional, Goedecker–Teter–Hutter pseudopotentials and TZV2P-GTH basis set. We applied the initial perturbation in the form of a $\delta$-pulse and we set the field strength parameter to 0.001 au. The RT-TDDFT time step was set to 0.00242 fs and we ran the simulation for up to 121 fs. We calculated the dipole moments of the whole simulation cell via the Berry phase approach for each RT-TDDFT step and we computed the polarizability tensors from the induced dipole moments. Using the [FFTW](https://www.fftw.org/) library, we applied fast fourier transformation to the polarizability tensors. Using the "plain-128" algorithm, we applied Pade approximants to the polarizabilities in the frequency domain.

Figure 5(a) shows the first absorption peak for naphthalene, generated from RT-TDDFT trajectories of different simulation lengths. RT-TDDFT trajectories with lengths of 121.0, 96.8, 72.6, 48.4, and 24.2 fs correspond to 50000, 40000, 30000, 20000, and 10000 RT-TDDFT steps, respectively. It is evident that, especially in the 24.2 fs trajectory, the absorption peak shifts to higher excitation energies due to insufficient data points and low resolution in the data set. The spectra seem to converge for simulation times greater than 100 fs. Figure 5(b) displays the absorption spectrum of naphthalene over the same excitation energy range, but this time with the inclusion of Padé approximants, extending the final number of data points to 80000 in each spectrum. The results indicate that, thanks to the increased resolution, the use of Padé approximants allows for a converged absorption spectrum even with RT-TDDFT simulation times as short as 24.2 fs, reducing the total computation time by approximately fivefold.

<div style="display:flex; justify-content: center; align-items: center;">
  <div style="width: 500px;">
  <img src="./img/absorption_naphthalene.png">
  <br>
  <div style="display: block; padding: 20px; color: gray; text-align: justify;">
   <b>Figure 5</b> Absorption spectra of the naphthalene molecule calculated from the RT-TDDFT trajectories of different simulation lengths a) without applying Padé approximants and b) applying Padé approximants using the "plain-128" algorithm.   
  </div>
  </div>
</div>

# Usage

There are two API functions that are needed in order to generate and evaluate a Thiele Padé interpolation. To create the Thiele Padé parameters call `create_thiele_Pade()` with the reference function arguments and values:
```fortran 
params_thiele = create_thiele_pade(n_par, x_ref, y_ref)
```
`x_ref`, `y_ref` must be of length `n_par`. After this step the parameters are stored in a fortran type called `params_thiele`. The parameters don't need to be accessed. In order to use the Padé model to evaluate function values with arbitrary function arguments you can use the API function `evaluate_thiele_Pade_at()`:
```fortran
y_return =  evaluate_thiele_pade_at(params_thiele, x_query)
```
`y_return` and `x_query` must be arrays of the same length. If the Padé model is not needed anymore, the parameters can be conveniently deallocated by:
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

### Available Options of Padé Interpolation
Fine-grained control over the generated pade model is provided by calling `create_thiele_pade()`with optional keyword arguments:
```fortran
params_thiele = create_thiele_pade(n_par, x_ref, y_ref,     &
                                   do_greedy = .true.,      &
                                   precision = 64,          &
                                   enforce_symmetry = "none")
y_return =  evaluate_thiele_pade_at(params_thiele, x_query)
```
The chosen options are stored in the model type and don't have to be repeated when the model is evaluated. All possible combinations of `do_greedy`, `precision` and `enforce_symmetry` options are supported. 

#### keyword argument `do_greedy` 
**Default:** `.true.` <br>
**Possible options:** `.true.`, `.false.`<br>
If true, a greedy algorithm is used to sort the reference points with the aim to lower the numerical noise. This comes at the cost of a slightly increased time to create the pade model. 

#### keyword argument `precision` 
**Default:** `128` if linked against GMP, else: `64`<br>
**Possible options:** any positive number greater zero of type `integer` <br>
The internal floating point precision in bit (not byte). Controls how floats are represented during creation and evaluation of the model using the GNU MP library for handling higher precision floats if a precision greater that 64 bit (double precision) is requested. The arrays containing the reference points (input) and also the evaluated function values (output) are in double precision independent of the `precision` keyword value. Note that a higher precision can increase the time of creating and evaluating the pade model drastically.

#### keyword argument `enforce_symmetry` 
**Default:** `none` <br>
**Possible options:** See table below. <br>
Force the Padé model to have a certain symmetry. If the symmetry of the underlying function is known, the user is advised to enforce this symmetry on the Padé model. This increases the predictive power of the model because more information about the function is provided.

| symmetry label | enforced symmetry | 
| --- | --- |
| `mirror_real` |  $f(z) = f(a+ib) = f(-a+ib)$ |
| `mirror_imag`  | $f(z) = f(a+ib) = f(a-ib)$ |
| `mirror_both` | $f(z) = f(a+ib) = f(a-ib) = f(-a+ib) = f(-a-ib) $ |
| `even` | $f(z) = f(-z)$ |
| `odd` | $f(z) = -f(-z)$ |
| `conjugate` | $f(z) = \overline{f(-z)}$ |
| `anti-conjugate` | $f(z) = -\overline{f(-z)}$ |
| `none` |  no symmetry enforced |


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
