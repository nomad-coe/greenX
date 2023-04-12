---
title: 'Time-frequency component of Green-X library: low scaling $\mathbf{GW}$ and $\mathbf{RPA}$ calculations'
tags:
  - FORTRAN
  - Low scaling $\mathbf{GW}$ calculations
  - Low scaling $\mathbf{RPA}$ calculations
  - Minimax approximation

authors:
  - name: Maryam Azizi
    orcid: 0000-0001-9089-1043
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Jan Wilhelm
    orcid: 0000-0003-0327-638X
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
date: 13 August 2017
bibliography: refs.bib
---

# Summary

Electronic structure calculations based on many-body perturbation theory (e.g. $\mathrm{GW}$ and random-phase approximation ($\mathrm{RPA}$)) involve manipulating functions
of time and frequency in the complex plane, like inhomogeneous Fourier transforms and analytic continuation from imaginary axis to real axis.
The time-frequency component of the Green-X library has the aim to provide such tools, in a coherent package.
In this initial release of the package, we provide a set of inhomogeneous time-frequency grids, as well as routines to compute the weights for calculating the susceptibility needed in low-scaling $\mathrm{RPA}$ and $\mathrm{GW}$ calculations, and for calculating the $\mathrm{GW}$ self-energy. They have been determined using the minimax approximation. Weights for both cosine and sine transformations are provided. 

The routines are freely available in GitHub, come with a build system, a set of tests, and related documentation. 

# Statement of need
$\mathrm{GW}$ calculations[@hedin1965new] have become part of the standard toolbox in computational condensed-matter physics for the calculation of photoemission and optoelectronic spectra of molecules and solids[@golze2019gw], [@reining2018gw], [@stankovski2011g], [@li2005quasiparticle], [@bruneval2008accurate]. There are, however, several algorithmic bottlenecks that render $\mathrm{GW}$ calculations challenging, especially for disordered systems with large simulations cells.

In conventional implementations of $\mathrm{GW}$ and $\mathrm{RPA}$, indeed, the computational cost of calculations increases with the fourth power of the system size $\mathrm{N}$, $\mathcal{O}(N^4)$. As a consequence, canonical $\mathrm{GW}$ calculations are usually limited to systems with a few hundred of atoms[@wilhelm2016gw], [@stuke2020atomic].
In order to reduce the computational cost of $\mathrm{GW}$ calculations, an alternative path relies on low-scaling algorithms that allow one to tackle larger and more realistic systems.
# References


