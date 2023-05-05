---
title: ' Time-frequency component of Green-X library: minimax grids for efficient RPA and GW calculations.'
tags:
  - FORTRAN
  - Low scaling GW calculations
  - Low scaling RPA calculations
  - Minimax approximation

authors:
  - name: Maryam Azizi
    orcid: 0000-0001-9089-1043
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Jan Wilhelm
    orcid: 0000-0001-8678-8246
    affiliation: 2
  - name: Dorothea Golze
    orcid: 0000-0002-2196-9350
    affiliation: 3
  - name: Matteo Giantomassi
    orcid: 0000-0002-7007-9813
    affiliation: 1
  - name: Ramón L. Panadés-Barrueta
    orcid: 0000-0003-4239-0978
    affiliation: 3
  - name: Francisco A. Delesma
    orcid: 0000-0001-6912-7745
    affiliation: 4
  - name: Alexander Buccheri
    orcid: 0000-0001-5983-8631
    affiliation: 5
  - name: Andris Gulans
    orcid: 0000-0001-7304-1952
    affiliation: 6
  - name: Patrick Rinke
    orcid: 0000-0003-1898-723X
    affiliation: 4
  - name: Claudia Draxl
    orcid: 0000-0003-3523-6657
    affiliation: 5
  - name: Xavier Gonze
    orcid: 0000-0002-8377-6829
    affiliation: 1
affiliations:
 - name: Université\ Catholique\ de\ Louvain,\ Louvain-la-Neuve, Belgium
   index: 1
 - name: University\ of Regensburg,\ Regensburg,\ Germany
   index: 2
 - name: Technische\ Universität\ Dresden,\ Dresden,\ Germany
   index: 3
 - name: Aalto University,\ Espoo,\ Finland
   index: 4
 - name: Humboldt-Universität\ zu\ Berlin,\ Berlin,\ Germany
   index: 5
 - name: University\ of\ Latvia,\ Riga,\ Latvia
   index: 6
date: 13 August 2017self
bibliography: refs.bib
---

# Summary

Electronic structure calculations based on many-body perturbation theory, such as the GW method and random-phase approximation (RPA), involve handling functions in the complex plane that vary over time and frequency. This includes performing inhomogeneous Fourier transforms. The Green-X library as a work package of the NOMAD CoE project aims to enable the electronic structure community to apply accurate beyond-DFT methodologies, based on Green's functions, to problems and systems currently out of reach with state-of-the-art supercomputers and software. The Green-X library's time-frequency component provides a coherent package of tools specifically designed for these tasks.

The initial release of the package comprises a collection of minimax time and frequency grids. The package also provides routines which compute the integration weights for the Fourier transform of the RPA susceptibility and the GW self-energy in low-scaling implementations. While the package targets low-scaling RPA and GW algorithms, its compact frequency grids can be used to reduce the computational prefactor in RPA implementations with conventional scaling. In addition, the time grids can be employed in Laplace-MP2 calculations.

The source code is freely available in GitHub, and comes equipped with a build system, a comprehensive set of tests, and detailed documentation.

# Statement of need

RPA is an accurate approach to compute the electronic correlation energy. The RPA correlation energy is non-local and includes long-range dispersion interactions [@ren2012random]. RPA takes dynamic electronic screening into account and is applicable to a wide range of systems, including molecules/clusters as well as extended structures [@eshuis2012electron],[@ren2012random]. The GW method [@hedin1965new] is based on the RPA susceptibility and has become the gold standard for the calculation of photoemission spectra of molecules and solids [@golze2019gw], [@reining2018gw], [@stankovski2011g], [@li2005quasiparticle], [@bruneval2008accurate], [@nabok2016accurate]. This can be followed by BSE (Bethe-Salpeter Equation) calculations of optoelectronic spectra [@onida2002electronic]. There are, however, several algorithmic bottlenecks that render RPA and GW calculations challenging, especially for disordered systems with large simulations cells.

In conventional implementations of RPA and GW, the computational cost increases with the fourth power of the system size $N$. As a consequence, canonical RPA and GW calculations are usually limited to systems with a few hundred of atoms [@delben2013electron], [@wilhelm2016gw], [@stuke2020atomic].
In order to reduce the computational cost of RPA and GW calculations, scaling reduction is a particularly promising strategy to tackle larger and more realistic systems.

Low-scaling algorithms rely on real space representation and time-frequency transforms, such as the real-space/imaginary-time approach [@rojas1995space]
with $\mathcal{O}(N^3)$ instead of $\mathcal{O}(N^4)$ complexity. Several cubic scaling GW algorithms have been recently implemented, e.g. in a plane-wave/projector-augmented-wave (PAW) GW code [@kaltak2014cubic], [@kaltak2014low], [@liu2016cubic] or with localized basis sets using Gaussian [@wilhelm2018toward], [@wilhelm2021low], [@duchemin2021cubic] or Slater-type orbitals [@forster2020low],[@foerster2021GW100], [@foerster2021loworder], [@foerster2023twocomponent]. Similarly, low-scaling RPA algorithms were implemented with different basis sets [@kaltak2014cubic], [@kaltak2014low], [@luenser2017], [@graf2018accurate], [@duchemin2019separable], [@drontschenko2022efficient].
An important property of low-scaling algorithms is the crossover point. The latter refers to the system size, where the low-scaling algorithm, which has usually a larger computational prefactor, becomes computationally more efficient than the canonical scheme [@wilhelm2018toward]. Another challenge for pioneering low-scaling GW algorithms was to reach
high numerical precision [@vlcek2017stochastic], [@wilhelm2018toward], [@forster2020low].

It has been found that the numerical precision of low-scaling GW algorithms is strongly coupled to the accuracy of the time-frequency grids as well as the technique used to numerically perform the analytic continuation of the self-energy matrix elements from the imaginary to the real axis. The existing implementations of the time-frequency and inverse transformations are concealed within electronic structure codes. Reusing these implementations elsewhere is either restricted by license requirements or dependencies on definitions made in the host code. In this work, we present an open-source package distributed under the Apache license (Version 2.0). The package handles the task by providing time and frequency grids and associated integration weights for low-scaling RPA and low-scaling GW as well as Laplace-Transform MP2 [@almloef1991elimination], [@jung2004scaled], [@kaltak2014low], [@glasbrenner2020efficient]. The package also offers an improvement over existing implementations by providing a broader selection of grids with an improved precision as shown in [@azizi_minimax]. As the generation of time-frequency grids can be numerically unstable and thus challenging and time demanding, the provided time-frequency grids could be considered as a comprehensive data set for low-scaling GW calculations for both molecules and solids.

# Mathematical framework

In the Adler-Wiser formula [@Adler1962],[@Wiser1963], the non-interacting susceptibility in the imaginary-frequency domain can be computed as follows

\begin{eqnarray}\label{susceptibility}
\chi^0( \bold{r}, \bold{r'}, i\omega ) = \sum_j^\text{occ}\sum_a^\text{unocc} \psi^*_a(\bold{r'})\psi_j(\bold{r'})\psi^*_j(\bold{r})\psi_a(\bold{r})\frac{2(\varepsilon_j - \varepsilon_a)}{\omega^2 + (\varepsilon_j - \varepsilon_a)^2},
\end{eqnarray}

where the indices $j$ and $a$ refer to occupied and unoccupied eigenstates, respectively, with energies $\varepsilon$ and wavefunctions $\psi$.

Writing Eq.(\ref{susceptibility}) in imaginary time 
\begin{eqnarray}\label{susceptibility_low}
    \hat{\chi}^0( \bold{r}, \bold{r'}, i\tau ) =\sum_j^\text{occ}\psi_j(\bold{r'})\psi^*_j(\bold{r})e^{\varepsilon_j|\tau|}\sum_a^\text{unocc}\psi^*_a(\bold{r'})\psi_a(\bold{r})e^{-\varepsilon_a|\tau|} 
\end{eqnarray}
allows to separate the two summations leading to a favorable computational scaling of $\mathcal{O}(N^3)$. 
Eq. \eqref{susceptibility_low} is employed as starting point of the low-scaling GW space-time method [@rojas1995space]. The imaginary-frequency result $\chi^0( \bold{r}, \bold{r'}, i\omega )$ is needed in RPA and GW calculations, thus requiring the Fourier transform
\begin{eqnarray}
    \chi^0( \bold{r}, \bold{r'}, i\omega ) = \int\limits_{-\infty}^{\infty} e^{-i \omega \tau} \, \hat{\chi}^0( \bold{r}, \bold{r'}, i\tau ) \, d\tau\label{ft_chi} \,.
\end{eqnarray}
The time integration is performed numerically using a discrete time and frequency grid [@rojas1995space], [@kaltak2014low]. For the sake of computational efficiency, the number of time and frequency grid points should be as small as possible. 
However, the function $\hat{\chi}^0( \bold{r}, \bold{r'}, i\tau )$ usually has long tails as well as very localized features in imaginary time. As such, the usual Fast Fourier Transform, with homogeneously spaced points, applied to Eq. \eqref{ft_chi} needs many sampling points. Instead, a nonuniform Fourier transform approach allows to significantly reduce the number of time and frequency grid points.

As a next step in GW, the screened Coulomb interaction $W( \bold{r}, \bold{r'}, i\omega )$ is computed which is then transformed to imaginary time
\begin{eqnarray} \label{ft_W}
       \hat{W}( \bold{r}, \bold{r'}, i\tau ) = \frac{1}{2\pi}\int\limits_{-\infty}^{\infty} e^{i \omega \tau} \, W( \bold{r}, \bold{r'}, i\omega )\, d\omega\,.
\end{eqnarray}
Moreover, a Fourier transform of the self-energy from imaginary time to imaginary frequency needs to be performed in the GW space-time method [@rojas1995space], [@liu2016cubic].

For constructing the discrete time and frequency grids, we split the computation of the functions $\hat{F}(i\tau)$ and $F(i\omega)$ in an even and an odd part [@liu2016cubic]
\begin{align}
\hat{F}(i\tau) &= \hat{F}^\text{even}(i\tau) + \hat{F}^\text{odd}(i\tau)
\\
{F}(i\omega) &=  {F}^\text{even}(i\omega) +  {F}^\text{odd}(i\omega)
\end{align}
with
\begin{align}
\hat{F}^\text{even}(i\tau)=\hat{F}^\text{even}(-i\tau) \hspace{2.4em} &\text{and} \hspace{1.9em} F^\text{even}(i\omega)=F^\text{even}(-i\omega)\,, \\
\hat{F}^\text{odd}(i\tau)=-\hat{F}^\text{odd}(-i\tau) \hspace{2em} &\text{and} \hspace{2em} F^\text{odd}(i\omega)=-F^\text{odd}(-i\omega)\,.
\end{align}

The corresponding discretized Fourier transforms read [@liu2016cubic]
\begin{align}\label{ct_st_even}
    F^\text{even}(i\omega_k) &= \sum_{j=1}^{N} \delta_{kj} \mathrm{cos}(\omega_k\tau_j)\hat{F}^\text{even}(i\tau_j)
    \\
    \hat{F}^\text{even}(i\tau_j)& = \sum_{k=1}^{N} \eta_{jk} \mathrm{cos}(\tau_j\omega_k)F^\text{even}(i\omega_k)\,,
 \\
    F^\text{odd}(i\omega_k)& = i\sum_{j=1}^{N} \lambda_{kj} \mathrm{sin}(\omega_k\tau_j)\hat{F}^\text{odd}(i\tau_j)\,,\label{st_odd_t_to_w}
    \\
    \hat{F}^\text{odd}(i\tau_j)&= -i \sum_{k=1}^{N} \zeta_{jk} \mathrm{sin}(\tau_j\omega_k)F^\text{odd}(i\omega_k),\label{ct_st_odd}
\end{align}
where $N$ is the number of grid points, $\{\tau_j\}_{j=1}^N, \tau_j>0$ are time grid points, $\{\omega_k\}_{k=1}^N,\omega_k>0$ are frequency grid points and
$\{\delta_{kj}\}_{k,j=1}^N$, $\{\eta_{jk}\}_{k,j=1}^N$, $\{\lambda_{kj}\}_{k,j=1}^N$, $\{\zeta_{jk}\}_{k,j=1}^N$ are integration weights. 
The transform \eqref{ct_st_even} is needed in low-scaling RPA calculations [@kaltak2014low], while the transforms \eqref{ct_st_even}-\eqref{st_odd_t_to_w} are needed in low-scaling GW calculations [@liu2016cubic]. We note  that $\{\zeta_{jk}\}_{k,j=1}^N$ is neither required for low-scaling RPA nor GW calculations and thus is not returned by our package.  

It is possible to construct optimal time and frequency points and optimal integration weights assuming that the functions have a specific shape in time and frequency. For RPA and GW, the assumed shape is inspired by Eqs. \eqref{susceptibility} and \eqref{susceptibility_low}, [@kaltak2014low]
\begin{align}
\hat{F}(i\tau) &= \sum_{x\in [x_\text{min},x_\text{max}]} \alpha_x \,e^{-\tau x}\hspace{3.65em}\text{for }\tau > 0\,,
\\
F(i\omega) &= \sum_{x\in [x_\text{min},x_\text{max}]} \alpha_x \,\frac{2x}{\omega^2+x^2}\hspace{2em}\text{for }\omega > 0\,,
\end{align}
where $x_\text{min} = \text{min}(\varepsilon_{a}-\varepsilon_{j})$ is the energy gap ($a$ refers to an unoccupied state, and $j$ to an occupied state) and $x_\text{max} = \text{max}(\varepsilon_{a}-\varepsilon_{j})$ is the maximum eigenvalue difference. $\alpha_x$ are expansion coefficients. Minimax grids are determined such that the error of the discretized Fourier transforms \eqref{ct_st_even} - \eqref{ct_st_odd} is minimized. For further details, we refer to [@braess2012nonlinear], [@kaltak2014low], [@hackbusch2019computation], [@liu2016cubic], [@azizi_minimax]. A minimax grid consists of $\{\tau_j\}_{j=1}^N, \tau_j>0$,   $\{\omega_k\}_{k=1}^N,\omega_k>0$,  $\{\delta_{kj}\}_{k,j=1}^N$, $\{\eta_{jk}\}_{k,j=1}^N$, $\{\lambda_{kj}\}_{k,j=1}^N$, $\{\zeta_{jk}\}_{k,j=1}^N$ and is specific to the interval  $[x_\text{min}, x_\text{max}]$. It is convenient to tabulate minimax grids for an interval $[1,R]$ with $R=x_\text{max}/x_\text{min}$ because these minimax grids relate to the minimax grid of the interval $[x_\text{min}, x_\text{max}]$ by rescaling [@kaltak2014low], [@hackbusch2019computation]. Our library provides minimax grids for grid points $N$ ranging from 6 to 34 and for a wide range of values for $R$.

For the calculation of the RPA correlation energy, a frequency integral needs to be computed which is discretized using the frequency grid $\{\omega_{k}\}_{k=1}^N$ and integration weights $\{\gamma_k\}_{k=1}^N$ [@kaltak2014low], [@delben2015enabling]. For the calculation of the MP2 correlation energy, a time quadrature can be performed using the time grid $\{\tau_j\}_{j=1}^N$ and corresponding integration weights $\{\sigma_j\}_{j=1}^N$ [@kaltak2014low]. The weights $\{\gamma_k\}_{k=1}^N$ and $\{\sigma_j\}_{j=1}^N$ are also provided by our package.

# Structure of the library

The Green-X library [@GitHub], [@azizi_minimax] includes eight components each of which tackles a critical problem in electronic structure calculations. In this work, we focus on the `GX-common` and `GX-TimeFrequency` components. The former contains common functionality used by all remaining library components, such as error handling and unit conversion utilities. The latter includes the source code for the new time-frequency analysis functionality presented in this publication. It comprises an API directory, which provides the necessary utilities for accessing and utilizing the new functionality, a source directory with the implementation of the algorithms described on the previous section, and a test directory with scripts for verifying the correctness of the implementation. The relevant directory tree section is the following:

```bash
.
├── CMakeLists.txt
├── developers.md
├── Doxyfile
├── GX-common
│   ├── CMakeLists.txt
│   └── src
│       ├── constants.f90
│       ├── error_handling.f90
│       ├── kinds.f90
│       ├── lapack_interfaces.f90
│       └── unit_conversion.f90
├── GX-TimeFrequency
│   ├── api
│   │   ├── api_utilities.f90
│   │   └── gx_minimax.f90
│   ├── CITATION.cff
│   ├── CMakeLists.txt
│   ├── LICENSE-2.0.txt
│   ├── README.md
│   ├── src
│   │   ├── gx_common.h
│   │   ├── minimax_grids.F90
│   │   ├── minimax_omega.F90
│   │   ├── minimax_tau.F90
│   │   └── minimax_utils.F90
│   ├── test
│   │   ├── conftest.py
│   │   ├── test_gx_minimax_grid.f90
│   │   ├── test_gx_minimax_grid.py
│   │   └── test_gx_tabulate_minimax.py
│   └── utilities
│       └── gx_tabulate_minimax.F90
├── LICENSE-2.0.txt
└── README.md
```

The project is written mostly in the Fortran programming language, supporting versions of the standard from 2008 and above. Some additional files employed in testing and error handling are written in C and Python. Modern Fortran features are extensively utilized, including object-oriented programming and some of the latest intrinsic procedures that have been added to the standard of the language. We have developed a clear interface between the module and client code, promoting better modularity and reusability. Additionally, we use allocatable arrays and automatic deallocation to simplify the code and avoid memory-related issues, such as leaks and dangling pointers. The implementation is robust and reliable, as we use error handling techniques to signal and recover from exceptions.

The time and frequency grids $\{\tau_k\}_{k=1}^N$, $\{\omega_j\}_{j=1}^N$ and the integration weights $\{\gamma_k\}_{k=1}^N$, $\{\sigma_j\}_{j=1}^N$ have been precomputed and are hard coded in the module `minimax_tau.f90` and `minimax_omega.f90` for various ranges $R$ of interest for solids and molecules. The module `minimax_grids.f90` calculates the weights $\{\delta_{kj}\}_{k,j=1}^N$, $\{\eta_{jk}\}_{k,j=1}^N$, $\{\lambda_{kj}\}_{k,j=1}^N$ of cosine and sine transforms during the runtime of the program via $L^2$ minimization [@azizi_minimax]. `minimax_grids.f90` also provides the corresponding maximum error of the $L^2$ minimization. To evaluate the precision of a global forward cosine transformation followed by backward cosine transformations, in the module `minimax grids` we provide a measure of such error as $\Delta_\text{CT}=\max_{j'j} | \sum_{k} \eta_{j'k} \cos(\tau_{j'}\omega_k) \cdot \delta_{kj} \cos(\omega_k\tau_j) - (\mathbb{I})_{j'j}|$ with $\mathbb{I}$ being the identity matrix. `minimax_utils.f90` contains auxiliary procedures and data structures for the main minimax routines.

To ensure the reproducibility of our results and improve the development process, we have implemented regression tests using the pytest infrastructure, combined with CMake and CTest. The tests rigorously verify the accuracy of the frequency grids and integration weights, the forward-backward cosine error, and the maximum errors of the transformations. Further details on the library structure, contents, build system, and testing are available in the various README.md files at the top level of the [@GitHub] repository and the individual component directories. The paper [@azizi_minimax] will also provide results of accuracy tests and realistic benchmarks.

# Acknowledgements

This work is supported by the European Union's Horizon 2020 research and innovation program under the grant agreement N° 951786 (NOMAD CoE). J.W. and D.G. acknowledge funding by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) via the Emmy Noether Programme (Project No. 503985532 and 453275048, respectively).

# References
