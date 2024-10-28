---
title: 'Analytic continuation component of the GreenX library: robust Padé approximants with symmetry constrains'
tags:
  - FORTRAN
  - Analytic Continuation
  - Pade 
  - GW approximation

authors:
  - name: Moritz Leucke
    orcid: 0009-0003-4381-0935
    affiliation: 1
  - name: Ramón L. Panadés-Barrueta
    orcid: 0000-0003-4239-0978
    affiliation: 1
  - name: Ekin E. Bas
    orcid: 0000-0002-0110-4691
    affiliation: 1
  - name: Dorothea Golze
    orcid: 0000-0002-2196-9350
    affiliation: 1
affiliations:
 - name: Faculty of Chemistry and Food Chemistry, Technische Universität Dresden, 01062 Dresden, Germany
   index: 1
date: 28 October 2024
bibliography: refs.bib
---

# Summary

Analytic continuation extends the domain of a given complex-valued function to a broader region in the complex plane. This technique is commonly used in many fields, for example in quantum mechanical methods such as the $GW$ method or real-time propagation algorithms, to continue a function from the imaginary axis to the real axis. 

In this work, we present the analytic continuation component of the GreenX library (GX-AC), which provides a Fortran API for the use of Padé approximants with and without symmetry constrains. The component uses the Thiele Padé algorithm to create Padé approximants and uses multiple-precision floats in combination with a greedy algorithm to mitigate the numerical instabilities commonly associated with fitting Padé approximants. The GX-AC component is distributed under the Apache 2 license and freely available on GitHub.


# Statement of need

![Application of the GX-AnalyticContinuation component to a model function with two poles (top left), an RT-TDDFT UV-vis Absorption spectrum (top right), the $GW$ self energy (bottom left) and the $GW$ screened coulomb interaction (bottom right). More information about the functions that are presented here can be found on the [website of the GX-AC component](https://nomad-coe.github.io/greenX/gx_ac.html).](ac_overview.pdf)

Analytic continuation (AC) is used in various scientific fields where complex analysis is relevant, like mathematical function theory, engineering and theoretical physics/chemistry in e.g.  quantum mechanics [@golze2019gw], quantum field theory [@nekrasov2024analytic], numerical methods for solving differential equations [@lope2002analytic] and real-time propagation methods [@li2020real]. In the following, we discuss the four examples depicted in Figure 1. The first example, shown in the top left of Figure 1, involves the application of AC to model functions, which may include Gamma functions [@luke1975error], Zeta functions [@iriguchi2007estimation], and others. 

In quantum field theory, AC can be applied to the frequently arising, complex-valued Green's functions, like the Green's function of the Hubbard model [@schott2016analytic].
However, Green's functions also appear in ab-initio many-body perturbation theory methods like the $GW$ approximation.
The $GW$ method [@hedin1965new] is considered the method of choice for predicting band structures of solids as well as electron removal and addition energies of molecules, as measured in direct and indirect photoemission experiments [@golze2019gw]. The complex-valued self energy is a central quantity in the $GW$ method, computed as the convolution of the Green's function $G$ and the screened interaction $W$. AC is a frequently used tool for continuing the self energy from the imaginary to the real frequency axis in conventional scaling $GW$ implementations [@van2015gw; @ren2012resolution; @gonze2009abinit; @wilhelm2017periodic] and low-scaling implementations [@liu2016cubic; @wilhelm2018toward; @wilhelm2021low; @graml2024low; @forster2020low; @forster2021low; @forster2021gw100; @forster2023two]. More recently, AC has also been applied to the screened interaction [@cdwac; @duchemin2020robust; @friedrich2019tetrahedron; @voora2020molecular; @samal2022modeling; @springer1998first; @kehry2023robust; @duchemin2021cubic] to e.g. reduce the computational scaling associated with core-level excitations [@cdwac]. The AC of a self energy $\Sigma$ and a screened Coulomb interaction $W$ are depicted as the second and third example in the bottom panel of Figure 1.

Our fourth and final example is the usage of AC in real-time propagation algorithms, such as real-time time-dependent density functional theory (RT-TDDFT) [@li2020real]. RT-TDDFT yields, for example, access to the absorption spectra of molecules and solids via the complex-valued dynamic polarizability tensor. The resolution of the RT-TDDFT absorption spectrum depends on the simulation length. It has been shown that applying Padé approximants to the dynamic polarizability tensor is an effective strategy for achieving higher spectral resolution with much shorter simulation times [@bruner2016accelerated; @mattiat2018efficient].
An illustrative UV-vis absorption spectrum, with and without the use of AC, is shown in the top right of Figure 1.


AC of analytic (holomorphic) functions is typically performed by  approximating the function with a rational function in one domain of the complex plane, typically along the imaginary axis. According to the identity theorem, the resulting rational function can then be evaluated over a broader domain of the complex plane, for example, along the real axis. 
Padé approximants are an established choice for rational functions. Their flexibility enables the approximation of functions with complicated pole structures [@golze2019gw]. Pade approximants can be expressed by the ratio of two polynomials with arbitrary order, or alternatively by a continued fraction.


The GreenX library aims to provide a suite of common tools, such as AC, for electronic structure codes based on the $GW$ method. The previously published first component of the GreenX library is the TimeFrequency component [@azizi2023time]. It provides minimax time and frequency grids for Random Phase Approximation (RPA) and $GW$ methods that were validated in a comprehensive benchmark study [@azizi2024validation]. In this work we present the second component of the GreenX library, the GX-AnalyticContinuation (GX-AC) component, which has the Apache-2.0 license. It provides a Fortran API for analytic continuation using Padé rational functions that can be easily integrated into other Fortran projects. The component uses the Thiele reciprocal difference method [@ThielePade_original; @ThielePade_Milne; @ThielePade_Baker] to obtain the Padé coefficients. Although the primary focus of the GreenX library are $GW$-based methods, the GX-AC component is suitable for any application where AC with Padé approximants can be used. Extensive benchmarks and the full documentation of GX-AC component can be found on the [component's website](https://nomad-coe.github.io/greenX/gx_ac.html).

Generating Padé approximants is prone to numerical instabilities caused by rounding errors that are amplified in the numerous differences in the Thiele-Padé algorithm [@PadeInstable; @cuyt1988instability; @jones1974numerical; @Beach2000], we employ two strategies to address these numerical instabilities. The first approach is to use  multiple precision floating point arithmetic for the implementation of the Thiele algorithm, minimizing the numerical noise caused by rounding errors. We use the GNU Multiple Precision (GMP) library [@GMPlib] to handle the multiple-precision floats. The advantage of this library is that it provides highly optimized assembly code for most of the processors available. This approach allows us to exceed the 128-bit precision limit typically supported by standard Fortran compilers. The second strategy involves using a greedy algorithm for Thiele Padé approximants, that has been validated in previous work [@greedy_pade1; @greedy_pade2; @cdwac]. The greedy algorithm is used to rearrange the function arguments of the reference function in order to make the model numerically more stable.

Another feature of the GX-AC component is to force the Padé model to exhibit a certain symmetry. This ensures that the approximant has the same symmetry as the reference function in the case that the symmetry of the reference function is known in advance, e.g. the screened interaction in the $GW$ is an even function [@duchemin2020robust]. Additionally, the enforced symmetry helps to increase the quality of the Padé approximant because every point of a given reference function also accounts for symmetrical equivalent points. Even, odd, conjugate and anti-conjugate function symmetry is supported by the component at this point.

# State of the field
To the best of the authors' knowledge, there are no Padé AC implementations that provide symmetry constraints and floating point precision beyond 128 bit, as the GX-AC component does. An open-source project that provides a Fortran implementation is the Padé Approximants repository by Johan Schött [@JohanSchott]. However, it uses the Beach algorithm [@Beach2000] in a quadruple precision implementation to generate Padé approximants and it provides an executable binary, rather than a Fortran API. Additionally, several Python, Julia and R implementations exist that offer an API in the specified language [@Adler_Pade_Pade_Approximant_2015; @2020SciPy-NMeth; @Montmorency; @jjgoings; @bennosski; @mjp98], but they do not offer a Fortran or C API.


# Acknowledgements
The authors acknowledge funding by the Emmy Noether Program of the German Research Foundation (Project No. 453275048). The German Research Foundation is gratefully acknowledged also for the support within the CRC 1415 (Chemistry of Synthetic Two-Dimensional Materials, Project No. 417590517).

# References
