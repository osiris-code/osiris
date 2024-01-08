---
layout: single
classes: wide
title: Grid Diagnostics
permalink: /other/grid_diagnostics
usemathjax: true

sidebar:
  nav: "other"
---

Grid quantities in OSIRIS share a common interface that allows the user to output several different types of grid diagnostics for a given grid quantity. The types of general diagnostics available are the following:

- Full dump of grid quantity to disk
- Data reduction operation: spatial averaging, spatial envelope,
  lineouts and slices
- Time averaging

Time averaging can be combined with spatial reduction -- for example, the user may select to report a lineout of a time average of a given quantity.

## Types of Diagnostics

In each diagnostics section, the selection of grid diagnostics to output is defined by a list of strings in the form:

`reports = diag1, diag2, ...,`

Each of the diagnostics must be specified in the form:

`quant [, tavg] [, savg | , senv | line, xj, j1(, j2) | slice, xj, j1 ]`

where

- **quant** specifies the quantity to be saved. See each specific diagnostics section to find out which grid quantities are available for diagnostics.
- **tavg** reports a time average of the selected quantity. The number of time steps averaged is controlled by the *n_tavg* parameter described below. Time averages can be used together with the data reduction operations described below.
- **savg** reports a spatial average of the selected quantity. The number of cells to average in each direction are specified using the *n_ave* parameter described below. The original grid size doesn't need to be a multiple of *n_ave*, but the diagnostic will be faster in this case.
- **senv** reports a spatial envelope of the selected quantity, calculated by finding the maximum absolute value of the quantity in a given number of neighboring grid cells. The number of cells to sample in every direction is specified using the *n_ave* parameter described below. The original grid size doesn't need to be a multiple of *n_ave*, but the diagnostic will be faster in this case.
- **line** (not available in 1D) extracts a single lineout from the selected grid quantity. The parameters specifying the line direction and position must follow:
  - **xj** must be "x1", "x2", or "x3" specifying the direction of the line to be extracted.
  - **j1**, **j2** - (*j2* in 3D only) are integers specifying the grid cell coordinates (1-based indexing) crossed by the line, and refer to the coordinates other than the line direction.
- **slice** (3D only) extracts a 2D slice from the selected grid quantity. The parameters specifying the slice direction and position must follow:
  - **xj** must be "x1", "x2", or "x3" specifying the normal direction of the slice to be extracted.
  - **j1** is an integer specifying the grid cell coordinates (1-based indexing) of the slice in the normal direction.

The *savg*, *senv*, *line*, and *slice* options are mutually exclusive and only one can be used in each diagnostic. A detailed description of the spatial average/envelope diagnostics can be found in the [average and envelope grids](average_envelope)
section of the documentation.  Note that the **j1** and **j2** indices refer to integer indices in the global grid, which start at 1 (Fortran-based indexing).

## Frequency of Diagnostics

The frequency of data dumps is controlled by three parameters for both normal and time-averaged diagnostics:

- **ndump_fac** - controls the frequency of full grid diagnostics.
- **ndump_fac_ave** - controls the frequency of spatial average/envelope grid diagnostics.
- **ndump_fac_lineout** - controls the frequency of lineout/slice diagnostics.

These values are multiplied by the global `ndump` parameter defined in the `time_step` section to get the actual diagnostic frequency. For example, if `ndump = 2`, `ndump_fac = 10`, and `ndump_fac_lineout = 1`, then the full grid will be output every 20 timesteps and a lineout will be output every 2 timesteps.

## Time and spatial filter parameters

The time average, spatial average, and spatial envelope diagnostics are controlled by the following parameters:

- **n_ave** specifies the number of gridpoints in each direction to average over/sample for spatially averaged/envelope dumps.
- **n_tavg** specifies the number of time steps to be used when calculating the time-averaged diagnostics.

## Example

Here's an example of a diag_emf section using grid diagnostics for a 3D simulation:

```text
 diag_emf 
 {
   reports = "e3",                     ! full E3 field
             "e2, line, x1, 64, 96",   ! lineout of E2 field along x1 at postion ix2 = 64, ix3 = 96
             "e3, senv",               ! spatial envelope of E3 
             "b1, tavg",               ! time average of B1
             "b2, tavg, savg",         ! Spatial average of the time average of B2
             "b3, slice, x2, 32",      ! Slice of B3 along x1-x3, at position ix2 = 32
      
   ndump_fac = 20,                     ! do full grid diagnostics at every 20 timesteps (assuming ndump = 1): 
   ndump_fac_ave = 10,                 ! do average/envelope grid diagnostics at every 10 timesteps
   ndump_fac_lineout = 1,              ! do lineouts/slices at every timestep
    
   n_ave(1:3) = 2,2,2,                 ! average/envelope 8 cells (2x2x2)
   n_tavg = 5,                         ! average 5 iterations for time averaged diagnostics 
 }
```
