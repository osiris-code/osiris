# Grid Diagnostics

Starting with revision 357 diagnostics of grid quantities in OSIRIS
share a common interface that allows the user to do several different
types of grid diagnostics besides the simple dumping of the grid
quantity to disk. The types of diagnostics available are the following:

- Full dump of grid quantity to disk.
- Data reduction operation: spatial averaging, spatial envelope,
  lineouts and slices.
- Time averaging.

The time averaging diagnostic can be used together with data reduction:
the user may select to report a lineout of a time average of a given
quantity.

## Types of Diagnostics

In each diagnostics section, the selection of grid diagnostics to do is
defined by a list of strings in the form:

`reports = diag1, diag2, ...,`

Each of the diagnostics selected must be in the form:

`quant [, tavg] [, savg | , senv | line, xi, i1(, i2) | slice, xi, i1 ]`

Where:

- **quant** - specifies the quantity to be saved. See each diagnostics
  section to find out which quantities are available for diagnostics.
- **tavg** - "time average", reports a time average of the selected
  quantity. The number of time steps averaged is controlled by the
  *n_tavg* parameter described below. Please note that time averages can
  be used together with the data reduction operations described below,
  i.e. it is possible to report a lineout of a time average.
- **savg** - "spatial average", reports a spatial average of the
  selected quantity. The number of cells to average in every direction
  are specified using the *n_ave* parameter described below. The
  original grid size doesn't need to be a multiple of *n_ave*, but the
  diagnostic will be faster in this case.
- **senv** - "spatial envelope", reports a spatial envelope of the
  selected quantity, that is calculated by finding the maximum absolute
  value of the quantity in a number of neighboring grid cells. The
  number of cells to sample in every direction are specified using the
  *n_ave* parameter described below. The original grid size doesn't need
  to be a multiple of *n_ave*, but the diagnostic will be faster in this
  case.
- **line** - extracts a single line from the selected grid quantity (not
  available in 1D). The parameters specifying the line direction and
  position must follow:
  - **xi** - must be "x1", "x2", or "x3" specifying the direction of the
    line to be extracted.
  - **i1**, **i2** - (*i2* in 3D only) Are integers specifying the grid
    cell coordinates crossed by the line, and refer to the coordinates
    other than the line direction.
- **slice** - extracts a 2D slice from the selected grid quantity (only
  available in 3D). The parameters specifying the slice direction and
  position must follow:
  - **xi** - must be "x1", "x2", or "x3" specifying the normal to the
    slice be extracted.
  - **i1** - Is an integers specifying the grid cell coordinates of the
    slice, and refer to the coordinates of the normal direction.

The *savg*, *senv*, *line*, and *slice* options are mutually exclusive
and only one can be used in each diagnostic. A detailed description of
the spatial average / envelope diagnostics can be found in the
[average/envelope grids](:User_Guide:_Average/Envelope_Grids "wikilink")
section of the user guide.

## Frequency of Diagnostics

The frequency of data dumps is controlled by three parameters, that
control both normal and time averaged diagnostics:

- **ndump_fac** - controls the frequency of full grid diagnostics.
- **ndump_fac_ave** - controls the frequency of spatial average /
  envelope grid diagnostics.
- **ndump_fac_lineout** - controls the frequency of lineout / slice
  diagnostics.

As it is usual in OSIRIS these values are multiplied by the global
*ndump* parameter defined in the time_step section to get the actual
diagnostic frequency.

## Time and Spatial Average parameters

The time and spatial average, and also spatial envelopes are controlled
by the following:

- **n_ave** specifies the number of gridpoints on each direction to
  average over/sample for spatially averaged/envelope dumps.
- **n_tavg** specifies the number of time steps to be used when
  calculating the time averaged diagnostics specified.

## Example

Here's an example of a diag_emf section using grid diagnostics:

```text
 diag_emf 
 {
   reports = "e3",                     ! full E3 field
             "e2, line, x1, 64, 96",   ! lineout of E2 field along x1 at postion ix2 = 64, ix3 = 96
             "e3, senv",               ! spatial envelope of E3 
             "b1, tavg",               ! time average of B1
             "b2, tavg, savg",         ! Spatial average of the time average of B2
             "b3, slice, x2, 32",      ! Slice of B3 along x1-x3, at position ix2 = 32
      
   ndump_fac = 20,                     ! do full grid diagnostics at every 20 timesteps
   ndump_fac_ave = 10,                 ! do average/envelope grid diagnostics at every 10 timesteps
   ndump_fac_lineout = 1,              ! do lineouts/slices at every timestep
    
   n_ave(1:3) = 2,2,2,                 ! average/envelope 8 cells (2x2x2)
   n_tavg = 5,                         ! average 5 iterations for time averaged diagnostics 
 }
```

