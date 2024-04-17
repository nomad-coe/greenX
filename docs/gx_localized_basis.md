---
layout: default
title: Localized Basis Component
tagline: GreenX Localized Basis
description: The implementation of the separable resolution of the identity
---

## General
This library provides a list localized basis set procedure, commonly occurring in post-scf calculations like second-orded Moller-Plesset (MP2) perturbation theory, random phase approximation (RPA) and <em>G<sub>0</sub>W<sub>0</sub></em> approach.

At the current stage, localized basis component of the GreenX library provide the separable resolution of the identity (RI) fitting coefficients to be used for the calculation of exact exchange energy, canonical RPA and GW methodologies. Also it can be used for the cubic-scaling RPA and GW, this feature is under development.

## Structure of the library
Workflow of the separable RI and cubic-scaling RPA with separable RI implementations as part of the GX-LocalizedBasis component of GreenX.

<p align="center">
  <img src="./img/Localizedbasis_structure.png" alt="Localizedbasis_structure" width="700">
</p>

## Benchmark

We have benchmarked the separable RI approach in our [recent publication](https://doi.org/10.1063/5.0184406) and we found that the separable RI reproduces the global RI reference within 1.0 meV for atomization energies predicted from correlated method such as MP2, coupled cluster single and doubles (CCSD), RPA, renormalized second order perturbation theory (rPT2), and <em>G<sub>0</sub>W<sub>0</sub></em>. For the Thiel set, which contains small organic molecules, these results are summarized in the following figure:

<p align="center">
  <img src="./img/Localizedbasis_validation.png" alt="Localizedbasis_validation" width="400">
</p>
<p align="center">
<em> Mean absolute deviations [meV] of separable RI with respect to RI-V utilizing the Thiel benchmark set.
</em></p>

Further benchmark calculations on the [GW100 test set](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00453) and [disordered carbon cluster](https://doi.org/10.1021/acs.chemmater.1c04279) can be found in [DOI:10.1063/5.0184406](https://doi.org/10.1063/5.0184406)



<button onclick="goBack()">Go Back</button>

<script>
function goBack() {
  window.history.back();
}
</script>
