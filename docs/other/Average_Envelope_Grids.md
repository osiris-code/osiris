---
layout: single
classes: wide
title: Average / Envelope grids
permalink: /other/average_envelope
usemathjax: true

sidebar:
  nav: "other"
---

## Introduction

The Average/Envelope grid diagnostics allow the user to significantly reduce the amount of output data by reducing the number of grid points in the output file. When using these diagnostics the code will determine the value of 1 grid point of the output file using the value of the multiple grid points of the original data that overlap with that point. This value can be calculated in 2 ways:

* **Average grids**. The output value is a weighted average of the original grid points. The weight used corresponds to the part of the original cell that overlaps with the averaged cell. This allows for the use of an arbitrary number of averaging points; however, the code will be faster if the number of averaging points divides the original grid evenly.
* **Envelope grids**. The output value is the maximum absolute value of the original grid points that overlap with the envelope cells. These are not "perfect" envelopes, but oversampling and smoothing in post-processing will give very good results.

## Using

All global grids in OSIRIS, EMF, charge (per species and total) and electric current can be saved using this diagnostic. The procedure is similar for all of them. Here's an example for the electric field:

```text
diag_emf  {
  ndump_fac_ave = 20,
  n_ave(1:2) = 2,3,
  reports = "e1, savg" ,
}
```

This will enable diagnostics of the spatially averaged charge of the given species at every `20 * ndump` iterations. The code will average 2 cells in the x1 direction and 3 cells in the x2 direction. The exact grid size of the output file is determined by doing an integer division between the global grid size and the values specified in `n_ave`. If the original grid had 128 x 128, the averaged grid will have (128:2) x (128:3) = 64 x 42 points. The spatial limits of the averaged grid are the same as the original grid.

Please check the reference guide for details on using this diagnostic:

* [Grid Diagnostics](grid_diagnostics) - Grid diagnostics parameters
* [particles](../reference/particles) - Global charge diagnostics
* [diag_species](../reference/diag_species) - Species charge diagnostics
* [diag_current](../reference/diag_current)- Electrical current diagnostics
* [diag_emf](../reference/diag_emf) -  Electro-magnetic field diagnostics

## Limitations

The code uses a fixed-size buffer for all these calculations. The maximum size for this buffer is set at compile time in `os-dutil.f90` and it is currently set to 64 MB. The code will process this diagnostic in chunks of sizes up to this value, so there is no limit for the final file size. The code is limited however to splitting the last dimension of the diagnostic: this means that 1D diagnostics are unlimited, 2D diagnostics will have a maximum size for rows, and unlimited columns and 3D diagnostics will have a limited size for each x1x2 slice, and unlimited x3 size.

In practice, these limits are very far away (with 64MB in 3D each slice is limited to 16 million cells or 4096^2). If you do hit these limits either increase the number of averaging points or recompile the code increasing the diagnostics buffer maximum size.
