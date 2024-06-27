---
layout: default
title: GreenX Library
tagline: GreenX
description: Library for Many-body Greens Functions on HPC
---

# Overview

An open-source library that supports exascale implementations of Green's-function-based methodologies. Its layered design separates higher-level functionalities (distinguishing between code-independent and code-family-specific parts) from architecture-dependent numerical routines, common to all code families.

Currently, our library supports the following electronic structure methods:

- conventional and low-scaling RPA
- low-scaling \\(GW\\)
- Laplace-transformed direct MP2

# Components
- [Analytic Continuation](gx_ac.md)
- [Plane-wave component](gx_planewave.md)
- [LAPW component ](gx_lapw.md)
- [Numerical-orbital component ](gx_localized_basis.md)
- [q=0 Coulomb component](gx_q0.md)
- [Minimax Time-Frequency](gx_time_frequency.md)


# Technical Details
 GreenX is written in Fortran 2008. The functionality needed for testing and error handling is written in C and Python. In the following we provide a few more details for current and future developers of GreenX. 
- [Structure of the library](structure.md)
- [Coding conventions and regression testing](tests.md)
