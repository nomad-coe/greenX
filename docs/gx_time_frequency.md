---
layout: page
title: Time-Frequency component
tagline: GreenX Time-Frequency
description: Time-Frequency component
---

# Description

# Usage

# Benchmarks

## CH<sub>4</sub> RPA

In this test, we evaluate the RPA total energy of CH4 using a Gauss-Legendre grid, a modified Gauss-Legendre grid (so far the standard in FHI-aims and abinit), and minimax grids. An accuracy of 10^-6 eV is reached with 10 minimax grid points while the modified Gauss-Legendre grids requires 36 points for this accuracy.

[CH4 benchmark](./img/ch4_becnh.png)

Error differences of the total RPA energy [eV] of methane calculated using the Gauss-Legendre, modified Gauss-Legendre and minimax imaginary frequency grid points. These differences were calculated with respect to the lowest RPA energy obtained with 34 minimax grid points. The ground state energy was calculated using the PBE exchange correlation functional in combination of the Tier2 basis set. The global resolution of identity (RI-V) approach was used for the calculation of the exact exchange and RPA correlation energy. The auxiliary basis functions for the RI-V method were generated automatically on the fly.

## GW100

---
