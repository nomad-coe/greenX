---
layout: default
title: GreenX Library
tagline: GreenX
description: Library for Many-body Greens Functions on HPC
---

# Overview

A new open-source library that supports exascale implementations of Green's-function-based methodologies. Its layered design will separate higher-level functionalities (distinguishing between code-independent and code-family-specific parts) from architecture-dependent numerical routines, common to all code families.

Currently, our library supports the following electronic structure methods:

- conventional and low-scaling RPA
- low-scaling \\(GW\\)
- Laplace-transformed direct MP2

# Components
- [Minimax Time-Frequency](gx_time_frequency.md)
- [Analytic Continuation](gx_ac.md)

# Technical Details
 GreenX is written in Fortran 2008. Functionality needed for testing and error handling is written in C and Python. In the following we provide a few more details for current and future developers of GreenX. 
- [Structure of the library](structure.md)
- [Coding conventions and regression testing](tests.md)
