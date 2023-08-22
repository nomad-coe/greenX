---
layout: default
title: Structure of the library
tagline: Structure
description: Structure of the library
---

# Structure of GX-common and GX-TimeFrequency

In this section, we focus on the components of the library which are ready to use, i.e., the  GX-common and GX-TimeFrequency components. GX-common provides functionality for all library components, such as error handling and unit conversion utilities. GX-TimeFrequency provides an API directory for the time-frequency transformations, a source directory, and a test directory with scripts for verifying the implementation. The relevant directory tree section is:

```plaintext
  |- CMakeLists.txt
  |- developers.md
  |- Doxyfile
  |- GX-common
  |  |- CMakeLists.txt
  |  |- src
  |  |  |- constants.f90
  |  |  |- error_handling.f90
  |  |  |- kinds.f90
  |  |  |- lapack_interfaces.f90
  |  |  |- unit_conversion.f90
  |- GX-TimeFrequency
  |  |- api
  |  |  |- api_utilities.f90
  |  |  |- gx_minimax.f90
  |  |- CITATION.cff
  |  |- CMakeLists.txt
  |  |- LICENSE.txt
  |  |- README.md
  |  |- src
  |  |  |- gx_common.h
  |  |  |- minimax_grids.F90
  |  |  |- minimax_omega.F90
  |  |  |- minimax_tau.F90
  |  |  |- minimax_utils.F90
  |  |- test
  |  |  |- conftest.py
  |  |  |- test_gx_minimax_grid.f90
  |  |  |- test_gx_minimax_grid.py
  |  |  |- test_gx_tabulate_minimax.py
  |  |- utilities
  |  |  |- gx_tabulate_minimax.F90
  |- LICENSE.txt
  |- README.md
```

The time and frequency grids $$\{\tau_k(R)\}_{k=1}^n$$, $$\{\omega_j(R)\}_{j=1}^n$$ and the integration weights $$\{\gamma_k(R)\}_{k=1}^n$$, $$\{\sigma_j(R)\}_{j=1}^n$$ are tabulated (hard coded) in the module `minimax_tau.f90` and `minimax_omega.f90` for various ranges $$R$$ of interest for solids and molecules. The module `minimax_grids.f90` calculates the weights $$\{\delta_{kj}\}_{k,j=1}^n$$, $$\{\eta_{jk}\}_{k,j=1}^n$$, $$\{\lambda_{kj}\}_{k,j=1}^n$$ of cosine and sine transforms during the runtime of the program. 'minimax_grids.f90' also provides the duality error $$\Delta_{\text{CT}}$$. `minimax_utils.f90` contains auxiliary procedures and data structures for the main minimax routines.

<button onclick="goBack()">Go Back</button>

<script>
function goBack() {
  window.history.back();
}
</script>
