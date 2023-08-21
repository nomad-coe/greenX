---
title: ' Time-frequency component of the GreenX library: minimax grids for efficient RPA and \textit{GW} calculations.'
tags:
  - FORTRAN
  - Low scaling \textit{GW} calculations
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
 - name: Institute of Condensed Matter and Nanoscience, UCLouvain, B-1348 \ Louvain-la-Neuve,\ Belgium\newline
   index: 1
 - name: Institute of Theoretical Physics and Regensburg Center for Ultrafast Nanoscopy (RUN), University of Regensburg, D-93053 \ Regensburg,\ Germany\newline
   index: 2
 - name: Faculty of Chemistry and Food Chemistry, Technische Universität Dresden, 01062 \ Dresden,\ Germany\newline
   index: 3
 - name: Department of Applied Physics, Aalto University, P.O. Box 11100, 00076 \ Aalto,\ Finland\newline
   index: 4
 - name: Institut für Physik und Iris Adlershof, Humboldt-Universität zu Berlin, Zum Großen Windkanal 2, 12489 \ Berlin,\ Germany\newline
   index: 5
 - name: Department of Physics, University of Latvia, Jelgavas iela 3, Riga, LV-1004 Latvia
   index: 6
date: 22 May 2023self
bibliography: refs.bib
---

# Summary

The central objects in many-body electronic structure theory, such as the \textit{GW} method and the random-phase approximation (RPA), are defined in the complex frequency or time domain. We present here the GX-TimeFrequency component of the GreenX library, which provides grids and weights for imaginary time-frequency transformations needed for Green's function based objects. The GreenX library emerged from the NOMAD Center of Excellence, whose objective is to enable accurate Green's function based electronic structure theory calculations on state-of-the-art supercomputers.

The package comprises minimax time and frequency grids [@Takatsuka2008; @kaltak2014low; @liu2016cubic] and corresponding quadrature weights to numerically compute time and frequency integrals of the correlation energy as well as weights for Fourier transforms between time and frequency grids. While we target low-scaling RPA and \textit{GW} algorithms, its compact frequency grids can also be used to reduce the computational prefactor in RPA implementations with conventional scaling. In addition, the time grids can be employed in Laplace-transformed direct MP2 (LT-dMP2) calculations. The GreenX source code is freely available on GitHub, and comes equipped with a build system, a comprehensive set of tests, and detailed documentation.

# Statement of need

RPA is an accurate approach to compute the electronic correlation energy. It is non-local, includes long-range dispersion interactions and dynamic electronic screening, and is applicable to a wide range of systems from 0 to 3 dimensions [@eshuis2012electron; @ren2012random]. The \textit{GW} method [@hedin1965new] is based on the RPA susceptibility and has become the method of choice for the calculation of direct and inverse photoemission spectra of molecules and solids [@golze2019gw;@reining2018gw]. \textit{GW} forms the basis for Bethe-Salpeter Equation (BSE) calculations of optical spectra [@onida2002electronic], where the \textit{GW} results are used as input.

Despite their wide adoption, RPA and \textit{GW} face computational challenges, especially for large systems. Conventional RPA and \textit{GW} implementations scale with the fourth power of the system size $N$ and are therefore usually limited to systems of a few tens to at most a hundred atoms [@wilhelm2016gw; @stuke2020atomic]. To tackle larger and more realistic systems, scaling reductions present a promising strategy to decrease the computational cost. Such low-scaling algorithms utilize real-space representations and time-frequency transformations, such as the real-space/imaginary-time approach [@rojas1995space] that reduces the complexity to $\mathcal{O}(N^3)$. Several such cubic-scaling \textit{GW} algorithms have recently been implemented, e.g. in a plane-wave/projector-augmented-wave (PAW) \textit{GW} code [@kaltak2014low;@liu2016cubic,@kutepov2017linearized] or with localized basis sets using Gaussian [@wilhelm2018toward;@wilhelm2021low;@duchemin2021cubic;@Graml2023] or Slater-type orbitals [@forster2020low;@foerster2021GW100;@foerster2021loworder;@foerster2023twocomponent]. Similarly, low-scaling RPA algorithms were implemented with different basis sets [@kaltak2014cubic;@kaltak2014low;@wilhelm2016rpa;@graf2018accurate;@duchemin2019separable;@drontschenko2022efficient].

An important consideration for low-scaling algorithms is the crossover point. Due to their larger pre-factor, low-scaling algorithms are typically more expensive for smaller systems and only become more cost-effective than canonical implementations for larger systems due to their reduced scaling [@wilhelm2018toward]. Furthermore, the numerical precision of low-scaling \textit{GW} algorithms is strongly coupled to the time-frequency treatment [@wilhelm2021low]. Early low-scaling \textit{GW} algorithms did not reach the same precision as canonical implementations [@vlcek2017stochastic; @wilhelm2018toward; @forster2020low]. Although appropriate Fourier transforms and corresponding time-frequency grids have been implemented [@liu2016cubic;@wilhelm2021low; @duchemin2021cubic; @foerster2021GW100], these implementations and grids are tied to particular codes and are often buried deeply inside the code. Furthermore, reuse of such implementations elsewhere is often restricted by license requirements or dependencies on definitions made in the host code.

In this work, we present the GX-TimeFrequency component of the GreenX library, an open-source package distributed under the Apache license (Version 2.0). GX-TimeFrequency provides time and frequency grids and the corresponding integration weights to compute correlation energies for Green's function implementations. It also provides Fourier weights to convert between imaginary time and imaginary frequency. The library can be used for low-scaling RPA and \textit{GW} implementations. BSE codes, which use (low-scaling) \textit{GW} as input, can also utilize our library. The minimax grids are also suitable for RPA implementations with conventional scaling [@delben2015enabling]: The minimax grids are more compact than, e.g., Gauss-Legendre grids, resulting in a reduction of the computational prefactor, while yielding the same accuracy [@delben2015enabling;@azizi_minimax]. We note that the minimax grids are not recommended for conventional imaginary-frequency-only \textit{GW} implementations [@ren2012resolution;@wilhelm2016gw] since they have not been optimized for the frequency integral of the self-energy.

While not being the main target of the library, the minimax time grids can also be utilized to calculate the LT-dMP2 correlation energy [@almloef1991elimination; @jung2004scaled; @Takatsuka2008; @kaltak2014low; @glasbrenner2020efficient]. The dMP2 term is one of two terms of the MP2 correlation energy. In a diagrammatic representation, dMP2 is the lowest order of the RPA correlation energy [@ren2012random]. The dMP2 correlation energy can be reformulated using the Laplace transform to obtain the LT-dMP2 expression which scales cubically in contrast to the $O(N^5)$ scaling of standard MP2.

# Mathematical framework

The single-particle Green's function $G$ and the non-interacting susceptibility $\chi^0$ are the starting point for a set of many-body perturbation theory (MBPT) methods. In canonical implementations, $\chi^0(\mathbf{r},\mathbf{r'},i\omega)$ is often expressed in the Adler-Wiser form [@Adler1962;@Wiser1963], where the sums over occupied (index $j$) and unoccupied (index $a$) single-particle states $\psi$ are coupled via their corresponding energies $\varepsilon$.

![Sketch of the methods supported by GX-TimeFrequency which start from $\hat{\chi}^0(i\tau)$. In addition to the discrete time and frequency grids $\{\tau_{j}\}$ and $\{\omega_{k}\}$, the library provides the corresponding weights $\{\sigma_{j}\}$ and $\{\gamma_{k}\}$ for the integration of the correlation energy $E_c$ as well as the Fourier weights $\delta_{kj}$, $\eta_{jk}$ and $\lambda_{kj}$ defined in Eqs. \eqref{ct_st_even}-\eqref{st_odd_t_to_w}. The bare and screened Coulomb interactions are indicated by $v(\mathbf{r},\mathbf{r}')=1/|\mathbf{r}-\mathbf{r}|'$ and $W(i\omega)$, respectively. $\epsilon(i\omega)$ is the dynamical dielectric function, $\Sigma$ the \textit{GW} self-energy, and AC stands for analytic continuation.\label{fig:flowchart}](flowchart.png)

The Adler-Wiser expression of $\chi^0(i\omega)$ can be transformed into the imaginary time domain, $\hat{\chi}^0(\mathbf{r},\mathbf{r'},i\tau)=-i G(\mathbf{r},\mathbf{r'},i\tau)G(\mathbf{r'},\mathbf{r},-i\tau)$, yielding the equation in the yellow box in Figure \ref{fig:flowchart}, where the two sums are separated leading to a favorable $\mathcal{O}(N^3)$ scaling. The polarizability $\hat{\chi}^0(i\tau)$ is the starting point for LT-dMP2 and low-scaling RPA and \textit{GW}, as summarized in Figure \ref{fig:flowchart}. The low-scaling \textit{GW} procedure shown in Figure \ref{fig:flowchart} is known as the space-time method and is given here in its original formulation for plane-wave codes [@rojas1995space].

The time-frequency integrals in Figure \ref{fig:flowchart} are performed numerically. All three methods in Figure \ref{fig:flowchart} require a discrete time grid $\{\tau_j\}_{j=1}^n$, where $n$ is the number of grid points. RPA and \textit{GW} additionally need the discrete frequency grid $\{\omega_{k}\}_{k=1}^n$. Since $\hat{\chi}^0( \mathbf{r}, \mathbf{r'}, i\tau )$ is sharply peaked around the origin and then decays slowly, homogeneous time and frequency grids are inefficient. For this reason, non-uniform grids like Gauss-Legendre [@rieger1999gw], modified Gauss-Legendre [@ren2012resolution] and the here presented minimax [@kaltak2014low] grids are used. The minimax grids include also integration weights for the computation of the correlation energies. For the calculation of the LT-dMP2 correlation energy $E_c^{\text{dMP2}}$ [@Takatsuka2008;@kaltak2014low], a time quadrature is performed, for which our library provides the integration weights $\{\sigma_j\}_{j=1}^n$. Similarly, the RPA correlation energy $E_c^{\text{RPA}}$ [@kaltak2014low;@delben2015enabling] is computed from a frequency quadrature using the integration weights $\{\gamma_k\}_{k=1}^n$.

Low-scaling RPA and \textit{GW} algorithms include the Fourier transform of $\hat{\chi}^0(i\tau)$ to $\chi^0(i\omega)$ (blue dashed box in Figure \ref{fig:flowchart}). The \textit{GW} space-time method performs two additional Fourier transforms: The screened Coulomb interaction $W(i\omega)$ is transformed to imaginary time (red dashed box), and the self-energy $\widehat{\Sigma}(i\tau)$ is Fourier transformed back to the imaginary frequency domain (green dashed box).

To convert between imaginary time and frequency grids, we introduce a nonuniform discrete cosine and sine transformation for even and odd functions $F^\text{even/odd}$, respectively [@liu2016cubic]. If the function $F$ is neither odd nor even we split the computation of functions $\hat{F}(i\tau)$ and $F(i\omega)$ into even and an odd parts [@liu2016cubic]

\begin{align}
\hat{F}(i\tau) = \hat{F}^\text{even}(i\tau) + \hat{F}^\text{odd}(i\tau) \hspace{1.9em} &\text{and} \hspace{1.9em} 
{F}(i\omega) =  {F}^\text{even}(i\omega) +  {F}^\text{odd}(i\omega)\,,\label{Fw_split}
\end{align}
with
\begin{align}
\hat{F}^\text{even}(i\tau)=\hat{F}^\text{even}(-i\tau) \hspace{2.4em} &\text{and} \hspace{1.9em} F^\text{even}(i\omega)=F^\text{even}(-i\omega)\,, \\
\hat{F}^\text{odd}(i\tau)=-\hat{F}^\text{odd}(-i\tau) \hspace{2em} &\text{and} \hspace{2em} F^\text{odd}(i\omega)=-F^\text{odd}(-i\omega)\,.
\end{align}

The corresponding discrete Fourier transforms read [@liu2016cubic]
\begin{align}\label{ct_st_even}
    F^\text{even}(i\omega_k) &= \sum_{j=1}^{n} \delta_{kj} \mathrm{cos}(\omega_k\tau_j)\hat{F}^\text{even}(i\tau_j)\,,
    \\
    \hat{F}^\text{even}(i\tau_j)& = \sum_{k=1}^{n} \eta_{jk} \mathrm{cos}(\tau_j\omega_k)F^\text{even}(i\omega_k)\,,\label{ct_even_w_to_t}
 \\
    F^\text{odd}(i\omega_k)& = i\sum_{j=1}^{n} \lambda_{kj} \mathrm{sin}(\omega_k\tau_j)\hat{F}^\text{odd}(i\tau_j)\,,\label{st_odd_t_to_w}
    \\
    \hat{F}^\text{odd}(i\tau_j)&= -i \sum_{k=1}^{n} \zeta_{jk} \mathrm{sin}(\tau_j\omega_k)F^\text{odd}(i\omega_k),\label{ct_st_odd}
\end{align}
where $\{\tau_j\}_{j=1}^n, \tau_j\,{>}\,0$ are again the time grid points, $\{\omega_k\}_{k=1}^n,\omega_k\,{>}\,0$ frequency grid points and $\{\delta_{kj}\}_{k,j=1}^n$, $\{\eta_{jk}\}_{k,j=1}^n$,$\{\lambda_{kj}\}_{k,j=1}^n$, $\{\zeta_{jk}\}_{k,j=1}^n$ the corresponding Fourier integration weights. $\hat{\chi}^0(i\tau)$ is an even function and we need thus the transform defined in Eq. \eqref{ct_st_even} to obtain $\chi^0(i\omega)$. The screened Coulomb interaction is also even and we use expression \eqref{ct_even_w_to_t} to convert $W(i\omega)$ to $\widehat{W}(i\tau)$. The self-energy is neither odd nor even and we use Eq. \eqref{Fw_split} in combination with Eqs. \eqref{ct_st_even} and  \eqref{st_odd_t_to_w} to transform $\widehat{\Sigma}(i\tau)$ to $\Sigma(i\omega)$ [@liu2016cubic]. The transformation defined in Eq. \eqref{ct_st_odd} is not required for the methods summarized in Figure \ref{fig:flowchart}, but is added for the sake of completeness and clarity.

Ideal grid parameters $\tau_j$, $\sigma_j$, $\omega_k$, $\gamma_k$, $\delta_{kj}$, $\eta_{jk}$, $\lambda_{kj}$ feature a vanishing error for the LT-dMP2 and RPA correlation energy integrations and Fourier transforms of $\chi^0,W$ and $\Sigma$ (Figure \ref{fig:flowchart}). We compute minimax grid parameters $\tau_j$, $\sigma_j$, $\omega_k$, $\gamma_k$ that minimize the maximum error of the LT-dMP2 and RPA correlation energy integration (Fig. \ref{fig:flowchart}) over all possible functions $\hat{\chi}^0( \mathbf{r}, \mathbf{r'}, i\tau )$ and $\chi^0( \mathbf{r}, \mathbf{r'}, i\omega )$
[@Takatsuka2008;@kaltak2014low;@liu2016cubic]. For this minimax grid optimization, we use the Remez algorithm [@kaltak2014low] which is an iterative, numerically ill-conditioned procedure requiring high numerical precision. As the generation of the minimax parameters $\tau_j$, $\sigma_j$, $\omega_k$, $\gamma_k$ is tedious, the most feasible strategy is to tabulate the computed minimax parameters $\{\tau_j\}_{j=1}^n$, $\{\sigma_j\}_{j=1}^n$, $\{\omega_k\}_{k=1}^n$, $\{\gamma_k\}_{k=1}^n$ for their later use in LT-dMP2, RPA, and \textit{GW} calculations.

It has been shown that minimax time and frequency grids $\{\tau_j\}_{j=1}^n$, $\{\omega_k\}_{k=1}^n$ are also suitable for performing Fourier transforms of $\chi^0,W$ and $\Sigma$ [@liu2016cubic]. With the knowledge of the tabulated $\{\tau_j\}_{j=1}^n$ and $\{\omega_k\}_{k=1}^n$ parameters, it is possible to use least-squares optimization to calculate the Fourier integration weights $\delta_{kj}$, $\eta_{jk}$, $\lambda_{kj}$ [@azizi_minimax]. The least-squares optimization can be executed by simple non-iterative linear matrix algebra which is straightforward and is done during the run time of the GreenX library [@azizi_minimax].

The optimal grid parameters $\tau_j$, $\sigma_j$, $\omega_k$,
$\gamma_k$, $\delta_{kj}$, $\eta_{jk}$, $\lambda_{kj}$ depend on the energy gap $\text{min}(\varepsilon_{a}-\varepsilon_{j})$ and the maximum eigenvalue difference $\text{max}(\varepsilon_{a}-\varepsilon_{j})$ of the material. We generated minimax grid parameters $\tau_j$, $\sigma_j$, $\omega_k$, $\gamma_k$ assuming energy differences $\varepsilon_a-\varepsilon_j\in [1,R]$, see details in Refs. [@kaltak2014low;@hackbusch2019computation;@azizi_minimax]. Our library stores minimax grid parameters $\{\tau_j(R)\}_{j=1}^n$, $\{\sigma_j(R)\}_{j=1}^n$, $\{\omega_k(R)\}_{k=1}^n$, $\{\gamma_k(R)\}_{k=1}^n$ for $n$ ranging from 6 to 34 and for different values of the *range* $R$ (on average 15 $R$-values for each $n$). For a material with energy gap $\Delta_\text{min}:=\text{min}(\varepsilon_{a}-\varepsilon_{j})$ and maximum eigenvalue difference $\Delta_\text{max}:=\text{max}(\varepsilon_{a}-\varepsilon_{j})$, one easily obtains the material-targeted minimax parameters $\{\tau_j^\text{mat}\}_{j=1}^n$, $\{\sigma_j^\text{mat}\}_{j=1}^n$, $\{\omega_k^\text{mat}\}_{k=1}^n$, $\{\gamma_k^\text{mat}\}_{k=1}^n$ from rescaling stored parameters with a range $R\ge\Delta_\text{max}/\Delta_\text{min}$ [@kaltak2014low;@hackbusch2019computation],
\begin{align}
  \omega_k^\text{mat} = \Delta_\text{min}\,\omega_k(R)\,,
 \hspace{1em}
 \gamma_k^\text{mat} = \Delta_\text{min}\,\gamma_k(R)\,,
 \hspace{1em}
\tau_j^\text{mat} = \frac{\tau_j(R)}{2\Delta_\text{min}}\,,
 \hspace{1em}
\sigma_j^\text{mat} = \frac{\sigma_j(R)}{2\Delta_\text{min}}\,.\label{eq:rescaling}
\end{align}

# Required input and output

GX-TimeFrequency requires as input the grid size $n$, the minimal eigenvalue difference $\Delta_{\text{min}}$, and the maximal eigenvalue difference $\Delta_{\text{max}}$. The output parameters are summarized in Table \ref{tab:output}. The library component retrieves the tabulated minimax parameters $\{\tau_j(R)\}_{j=1}^n$, $\{\sigma_j(R)\}_{j=1}^n$, $\{\omega_k(R)\}_{k=1}^n$, $\{\gamma_k(R)\}_{k=1}^n$ of the requested grid size $n$ for the smallest range $R$ that satisfies $R \ge \Delta_\text{max}/\Delta_\text{min}$. GX-TimeFrequency then rescales the retrieved minimax parameters according to Equation \eqref{eq:rescaling} with $\Delta_\text{min}$ and prints the results $\{\tau_j^\text{mat}\}_{j=1}^n$, $\{\sigma_j^\text{mat}\}_{j=1}^n$, $\{\omega_k^\text{mat}\}_{k=1}^n$, $\{\gamma_k^\text{mat}\}_{k=1}^n$. Fourier integration weights are computed on-the-fly via least-squares optimization. To evaluate the precision of a global forward cosine transformation followed by backward cosine transformations, we provide a measure of such error as 
\begin{align}
\Delta_\text{CT}=\max_{j,j'\in\{1,2,\ldots,n\}} \left| \sum_{k=1}^n \eta_{j'k} \cos(\tau_{j'}\omega_k) \cdot \delta_{kj} \cos(\omega_k\tau_j) - (\mathbb{I})_{j'j}\right|
\end{align} 
with $\mathbb{I}$ being the identity matrix. Inputs and outputs are in atomic units.

| Output| Description|Methods using the output|Computation|
|---|--------|---------|-------|
|$\{\tau_j^\text{mat}\}_{j=1}^n$&nbsp; &nbsp;  | time points | LT-dMP2, ls RPA, ls \textit{GW} | tabulated + rescaling |
|$\{\sigma_j^\text{mat}\}_{j=1}^n$ | time integration weights | LT-dMP2  | tabulated + rescaling
|$\{\omega_k^\text{mat}\}_{k=1}^n$ | frequency points | ls & canonical RPA, ls \textit{GW}  | tabulated + rescaling |
|$\{\gamma_k^\text{mat}\}_{k=1}^n$ | freq. integration weights | ls & canonical RPA | tabulated + rescaling  |
|$\{\delta_{kj}\}_{k,j=1}^n$ | Fourier weights | ls RPA, ls \textit{GW} | on-the-fly L2 opt  
|  $\{\eta_{jk}\}_{k,j=1}^n$ | Fourier weights | ls \textit{GW} | on-the-fly L2 opt|
|$\{\lambda_{kj}\}_{k,j=1}^n$ | Fourier weights | ls \textit{GW} | on-the-fly L2 opt  
|  $\Delta_\text{CT}$ | duality error cosine transforms | ls \textit{GW} | on-the-fly  |
: Output returned by the GX-TimeFrequency component of GreenX. We abbreviate low-scaling as ls, and least-squares optimization as L2 opt.\label{tab:output}

# Acknowledgements

This work is supported by the European Union's Horizon 2020 research and innovation program under the grant agreement N° 951786 (NOMAD CoE). J.W. and D.G. acknowledge funding by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) via the Emmy Noether Programme (Project No. 503985532 and 453275048, respectively). A.G. acknowledges funding provided by the European Regional Development Fund via the Central Finance and Contracting Agency of Republic of Latvia (Project No. 1.1.1.5/21/A/004).

# References
