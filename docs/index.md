---
layout: default
title: GreenX Library
tagline: GreenX
description: Library for Many-body Greens Functions on HPC
---

# Overview

A new open-source library that supports exascale implementations of Green-function-based methodologies. Its layered design will separate higher-level functionalities (distinguishing between code-independent and code-family-specific parts) from architecture-dependent numerical routines, common to all code families. GreenX is written in Fortran 2008. Functionality needed for testing and error handling is written in C and Python. We utilize modern Fortran features such as object-oriented programming and intrinsic procedures that are available in Fortran 2008. We have developed a clear interface between the module code (our library) and the client code (MBPT code), promoting better modularity and reusability. Additionally, we use allocatable arrays and automatic deallocation to simplify the code and avoid memory-related issues, such as leaks and dangling pointers. The implementation is robust and reliable, as we use error handling techniques to highlight and recover from exceptions. 
- [Structure of the library](structure.md)
# Components
- [Minimax Time-Frequency](gx_time_frequency.md)
- [Analytic Continuation](gx_ac.md)
