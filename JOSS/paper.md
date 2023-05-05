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

The central objects in many-body electronic structure theory, such as the GW method and random-phase approximation (RPA), are defined in the complex frequency or time domain. Fourier transforms are used to transform between time an frequency. The Green-X library presented here provides the infrastructure and tools for such Fourier transforms with minimax grids. The Green-X library emerged from the NOMAD Center of Excellence, whose objective is to enable accurate Green's function based electronic structure theory calculations on  state-of-the-art supercomputers.

Green-X provides minimax time and frequency grids and the corresponding integration weights for Fourier transforms. While the package targets low-scaling RPA and GW algorithms, its compact frequency grids can also be used to reduce the computational prefactor in RPA implementations with conventional scaling. In addition, the time grids can be employed in Laplace-MP2 calculations. The Green-X source code is freely available on GitHub, and comes equipped with a build system, a comprehensive set of tests, and detailed documentation.

# Statement of need

RPA is an accurate approach to compute the electronic correlation energy. It is non-local, includes long-range dispersion interactions and dynamic electronic screening and is applicable to a wide range of systems from 0 to 3 dimensions [@eshuis2012electron],[@ren2012random]. The GW method [@hedin1965new] is based on the RPA susceptibility and has become the method of choice for the calculation of direct and inverse photoemission spectra of molecules and solids [@golze2019gw], [@reining2018gw], [@stankovski2011g], [@li2005quasiparticle], [@bruneval2008accurate], [@nabok2016accurate]. $GW$ forms the basis for the Bethe-Salpeter Equation (BSE) calculations of optical spectra [@onida2002electronic] and for advanced correlation methods like second-ordere screened exchange (SOSEX) [@Ren/etal:2015],[@Wang/Rinke/Ren:2021]


Despite their wide adoption, RPA, GW and BSE face computational challenges, especially for large systems. Conventional RPA and GW implementations scale with the fourth power of the system size $N$ and are therefore limited to systems of a few hundred atoms [@delben2013electron], [@wilhelm2016gw], [@stuke2020atomic]. 
To tackle larger and more realistic systems, scaling reductions present a promising strategy to reduce the computational cost. Such low-scaling algorithms utilize real space representation and time-frequency transformations, such as the real-space/imaginary-time approach [@rojas1995space] that recudes the complexity to $\mathcal{O}(N^3)$. Several such cubic scaling GW algorithms have recently been implemented, e.g. in a plane-wave/projector-augmented-wave (PAW) GW code [@kaltak2014cubic], [@kaltak2014low], [@liu2016cubic] or with localized basis sets using Gaussian [@wilhelm2018toward], [@wilhelm2021low], [@duchemin2021cubic] or Slater-type orbitals [@forster2020low],[@foerster2021GW100], [@foerster2021loworder], [@foerster2023twocomponent]. Similarly, low-scaling RPA algorithms were implemented with different basis sets [@kaltak2014cubic], [@kaltak2014low], [@luenser2017], [@graf2018accurate], [@duchemin2019separable], [@drontschenko2022efficient].
An important consideration for low-scaling algorithms is the crossover point. Due to their larger pre-factor, low scaling algorithms are typically more expensive for smaller systems and only become more cost effective than canonical implementations for larger systems due to their reduced scaling [@wilhelm2018toward]. Furthermore, early low-scaling GW algorithms did not reach the same numerical accuracy as canonical implementations [@vlcek2017stochastic], [@wilhelm2018toward], [@forster2020low]. The numerical precision of low-scaling GW algorithms is strongly coupled to the time-frequency treatment the analytic continuation of the self-energy matrix elements from the imaginary to the real frequency axis. Although appropriate Fourier transforms and corresponding time-frequency grids have been implemented, these implementations are tied to particular codes and are often burried deeply inside the code. Furthermore, reuse of such implementations elsewhere is often restricted by license requirements or dependencies on definitions made in the host code. 

In this work, we present the Green-X libary, an open-source package distributed under the Apache license (Version 2.0). Green-X provides time and frequency grids and corresponding integration weights for low-scaling Green's function implementations (e.g., RPA,GW, SOSEX, Laplace-Transform MP2)  [@almloef1991elimination], [@jung2004scaled], [@kaltak2014low], [@glasbrenner2020efficient]. As the generation of time-frequency grids can be numerically unstable and thus challenging and time demanding, the provided time-frequency grids could be considered as a comprehensive data set for low-scaling GW calculations for both molecules and solids.

# Mathematical framework

The central objects in many-body perturbation theory are the single-particle Green's function $G$, the non-interacting susceptibility $\chi^0$, the screened Coulomb interaction $W$ and the self-energy $\Sigma$. In canonical implementations, the non-interacting susceptibility is often expressed in the Adler-Wiser form [@Adler1962],[@Wiser1963]

\begin{eqnarray}\label{susceptibility}
\chi^0( \bold{r}, \bold{r'}, i\omega ) = \sum_j^\text{occ}\sum_a^\text{unocc} \psi^*_a(\bold{r'})\psi_j(\bold{r'})\psi^*_j(\bold{r})\psi_a(\bold{r})\frac{2(\varepsilon_j - \varepsilon_a)}{\omega^2 + (\varepsilon_j - \varepsilon_a)^2},
\end{eqnarray}

where the indices $j$ and $a$ refer to occupied and unoccupied single particle states $\psi$ with energies $\varepsilon$.

Transforming Eq.(\ref{susceptibility}) into the imaginary time domain 
\begin{eqnarray}\label{susceptibility_low}
    \hat{\chi}^0( \bold{r}, \bold{r'}, i\tau ) &=&\sum_j^\text{occ}\psi_j(\bold{r'})\psi^*_j(\bold{r})e^{\varepsilon_j|\tau|}\sum_a^\text{unocc}\psi^*_a(\bold{r'})\psi_a(\bold{r})e^{-\varepsilon_a|\tau|} \\
    &=& G(\mathbf{r},\mathbf{r'},i\tau)G(\mathbf{r'},\mathbf{r},-i\tau)
\end{eqnarray}
separates the two summations that can also be written as product of two Green's functions and leads to a favorable $\mathcal{O}(N^3)$ scaling.
Eq. \eqref{susceptibility_low} is typically the starting point for low-scaling GW implementations like the GW space-time method [@rojas1995space].

Once the non-interacting susceptibility has been calculated in imaginary time, it is Fourier transformed to imaginary frequency
\begin{eqnarray}
    \chi^0( \bold{r}, \bold{r'}, i\omega ) = \int\limits_{-\infty}^{\infty} e^{-i \omega \tau} \, \hat{\chi}^0( \bold{r}, \bold{r'}, i\tau ) \, d\tau\label{ft_chi} \,.
\end{eqnarray}
The time integration is performed numerically using a discrete time and frequency grid [@rojas1995space], [@kaltak2014low]. Since $\hat{\chi}^0( \bold{r}, \bold{r'}, i\tau )$ is sharply peaked around the origin and then decays slow, homogeneous time and frequency grids are inefficient. For this reason, non-uniform grids like Gauss-Legendre [@rieger1999gw], modified Gauss-Legendre [@ren2012resolution] and the here presented MiniMax grids are used. 

In GW, the screened Coulomb interaction $W( \bold{r}, \bold{r'}, i\omega )$ is then computed in imaginary frequency and subsequently transformed to imaginary time
\begin{eqnarray} \label{ft_W}
       \hat{W}( \bold{r}, \bold{r'}, i\tau ) = \frac{1}{2\pi}\int\limits_{-\infty}^{\infty} e^{i \omega \tau} \, W( \bold{r}, \bold{r'}, i\omega )\, d\omega\,.
\end{eqnarray}
Lastly, the self-energy is Fourier transformed back to imaginary frequency [@rojas1995space], [@liu2016cubic].

For the discrete time and frequency grids, we split the computation of functions $\hat{F}(i\tau)$ and $F(i\omega)$ into even and an odd parts [@liu2016cubic]
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
where $N$ is the number of grid points. $\{\tau_j\}_{j=1}^N, \tau_j>0$ are time grid points, $\{\omega_k\}_{k=1}^N,\omega_k>0$ frequency grid points and
$\{\delta_{kj}\}_{k,j=1}^N$, $\{\eta_{jk}\}_{k,j=1}^N$, $\{\lambda_{kj}\}_{k,j=1}^N$, $\{\zeta_{jk}\}_{k,j=1}^N$ the corresponding integration weights. 
The transform \eqref{ct_st_even} is needed in low-scaling RPA calculations [@kaltak2014low], while the transforms \eqref{ct_st_even}-\eqref{st_odd_t_to_w} are needed in low-scaling GW calculations [@liu2016cubic]. We note  that $\{\zeta_{jk}\}_{k,j=1}^N$ is neither required for low-scaling RPA nor GW calculations and thus is not returned by Green-X.  

It is possible to construct optimal time and frequency points and optimal integration weights assuming that the functions have a specific shape in time and frequency. For RPA and GW, the assumed shape is inspired by Eqs. \eqref{susceptibility} and \eqref{susceptibility_low}, [@kaltak2014low]
\begin{align}
\hat{F}(i\tau) &= \sum_{x\in [x_\text{min},x_\text{max}]} \alpha_x \,e^{-\tau x}\hspace{3.65em}\text{for }\tau > 0\,,
\\
F(i\omega) &= \sum_{x\in [x_\text{min},x_\text{max}]} \alpha_x \,\frac{2x}{\omega^2+x^2}\hspace{2em}\text{for }\omega > 0\,,
\end{align}
where $x_\text{min} = \text{min}(\varepsilon_{a}-\varepsilon_{j})$ is the energy gap ($a$ refers to an unoccupied state, and $j$ to an occupied state) and $x_\text{max} = \text{max}(\varepsilon_{a}-\varepsilon_{j})$ is the maximum eigenvalue difference. $\alpha_x$ are expansion coefficients. Minimax grids are determined such that the error of the discretized Fourier transforms \eqref{ct_st_even} - \eqref{ct_st_odd} is minimized. For further details, we refer to [@braess2012nonlinear], [@kaltak2014low], [@hackbusch2019computation], [@liu2016cubic], [@azizi_minimax]. A minimax grid consists of $\{\tau_j\}_{j=1}^N, \tau_j>0$,   $\{\omega_k\}_{k=1}^N,\omega_k>0$,  $\{\delta_{kj}\}_{k,j=1}^N$, $\{\eta_{jk}\}_{k,j=1}^N$, $\{\lambda_{kj}\}_{k,j=1}^N$, $\{\zeta_{jk}\}_{k,j=1}^N$ and is specific to the interval  $[x_\text{min}, x_\text{max}]$. It is convenient to tabulate minimax grids for an interval $[1,R]$ with $R=x_\text{max}/x_\text{min}$ because these minimax grids relate to the minimax grid of the interval $[x_\text{min}, x_\text{max}]$ by rescaling [@kaltak2014low], [@hackbusch2019computation]. Our library provides minimax grids for grid points $N$ ranging from 6 to 34 and for a wide range of values for $R$.

For the calculation of the RPA correlation energy, a frequency integral needs to be computed which is discretized using the frequency grid $\{\omega_{k}\}_{k=1}^N$ and integration weights $\{\gamma_k\}_{k=1}^N$ [@kaltak2014low], [@delben2015enabling]. For the calculation of the MP2 correlation energy, a time quadrature can be performed using the time grid $\{\tau_j\}_{j=1}^N$ and corresponding integration weights $\{\sigma_j\}_{j=1}^N$ [@kaltak2014low]. The weights $\{\gamma_k\}_{k=1}^N$ and $\{\sigma_j\}_{j=1}^N$ are also provided by our package.

# Structure of the library

The Green-X library [@GitHub], [@azizi_minimax] will eventually provide a variety of tools for advanced electronic structure calculations. In this work, we focus on the `GX-common` and `GX-TimeFrequency` components. `GX-common` provides functionality for all library components, such as error handling and unit conversion utilities. `GX-TimeFrequency` provides an API directory for the time-frequency transformations, a source directory, and a test directory with scripts for verifying the implementation. The relevant directory tree section is:

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

Green-X is written mostly in the Fortran programming language, supporting versions of the standard from 2008 and above. Functionality needed for testing and error handling is written in C and Python. We utilize modern Fortran features such as object-oriented programming and some of the latest intrinsic procedures that have been added to the standard of the language. We have developed a clear interface between the module and client code, promoting better modularity and reusability. Additionally, we use allocatable arrays and automatic deallocation to simplify the code and avoid memory-related issues, such as leaks and dangling pointers. The implementation is robust and reliable, as we use error handling techniques to signal and recover from exceptions.

The time and frequency grids $\{\tau_k\}_{k=1}^N$, $\{\omega_j\}_{j=1}^N$ and the integration weights $\{\gamma_k\}_{k=1}^N$, $\{\sigma_j\}_{j=1}^N$ have been precomputed and are hard coded in the module `minimax_tau.f90` and `minimax_omega.f90` for various ranges $R$ of interest for solids and molecules. The module `minimax_grids.f90` calculates the weights $\{\delta_{kj}\}_{k,j=1}^N$, $\{\eta_{jk}\}_{k,j=1}^N$, $\{\lambda_{kj}\}_{k,j=1}^N$ of cosine and sine transforms during the runtime of the program via $L^2$ minimization [@azizi_minimax]. `minimax_grids.f90` also provides the corresponding maximum error of the $L^2$ minimization. To evaluate the precision of a global forward cosine transformation followed by backward cosine transformations, in the module `minimax grids` we provide a measure of such error as $\displaystyle \Delta_\text{CT}=\max_{j'j} | \sum_{k} \eta_{j'k} \cos(\tau_{j'}\omega_k) \cdot \delta_{kj} \cos(\omega_k\tau_j) - (\mathbb{I})_{j'j}|$ with $\mathbb{I}$ being the identity matrix. `minimax_utils.f90` contains auxiliary procedures and data structures for the main minimax routines.

To ensure the reproducibility of our results and improve the development process, we have implemented regression tests using the pytest infrastructure, combined with CMake and CTest. The regression tests rigorously verify the accuracy of the frequency grids and integration weights, the forward-backward cosine error, and the maximum errors of the transformations. Further details on the library structure, contents, build system, and testing are available in the various README.md files at the top level of the [@GitHub] repository and the individual component directories. Further benchmarks will be provided elsewhere [@azizi_minimax].

# Acknowledgements

This work is supported by the European Union's Horizon 2020 research and innovation program under the grant agreement N° 951786 (NOMAD CoE). J.W. and D.G. acknowledge funding by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) via the Emmy Noether Programme (Project No. 503985532 and 453275048, respectively).

# References
