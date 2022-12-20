---
layout: single
classes: wide
title: Electric Current Diagnostics
permalink: /reference/diag_current
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the electric current diagnostic settings, used for reporting the total simulation current -- the sum of currents from all species.

The current for individual species can separately be output as reports in the [`diag_species`](Species_Diagnostics.md) namelist. There are two other subtle differences between the individual current diagnostics and the total current diagnostic: 1) diagnostics for the total current include optionally-applied current smoothing and solver-dependent current corrections (see [`emf_solver`](Electro-Magnetic_Field_Solver.md)), while these are absent for the individual species diagnostics; and 2) currents for individual species are not deposited using the same current discretization scheme as the total current. For these reasons, diagnostics for the current of individual species output by `diag_species` do not exactly add to the total current output by `diag_current`; however, for almost all problems, the individual species currents are perfectly valid for analysis.

This section is optional; if not present, the code will not do any electrical current diagnostics. It accepts the following parameters:

- **ndump_fac**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **prec**, integer, default = 4
- **n_tavg**, integer, default = -1
- **n_ave**(x_dim), integer, default = -1
- **reports**(:), character(\*), default = "-"

**ndump_fac** controls the frequency of full grid diagnostics. This value is multiplied by `ndump` value specified in the `time_step` section to determine the number of iterations between each diagnostic dump. If set to 0 the writing of this diagnostic information is disabled.

**ndump_fac_ave** controls the frequency of spatial average / envelope grid diagnostics. This value is multiplied by the `ndump` value specified in the `time_step` section to determine the number of iterations between each diagnostic dump. If set to 0 the writing of this diagnostic information is disabled.

**ndump_fac_lineout** controls the frequency of lineout / slice diagnostics. This value is multiplied by the `ndump` value specified in the `time_step` section to determine the number of iterations between each diagnostic dump. If set to 0 the writing of this diagnostic information is disabled.

**n_tavg** specifies the number of time steps to be used when calculating the time averaged diagnostics. The frequency of these diagnostics is controlled by the `ndump_fac` parameter described above.

**n_ave** number of gridpoints on each direction to average over for spatially averaged dumps. The frequency of these diagnostics is controlled by the `ndump_fac_ave` parameter described above.

**prec** controls the numerical precision used for grid diagnostics. The default is to save data in single precision (prec = 4). If the user wants data to be saved in double precision this parameter must be set to 8. This option is ignored if OSIRIS is compiled in single precision.

**reports** specifies the grid quantities to report, including spatial/time averaging, lineouts, etc., as described in the [grid diagnostics section](:Reference_Guide:_Grid_Diagnostics "wikilink"). The available quantities are:

- "j1", "j2", "j3" - total current components
- "div_j" - divergence of total current

Here is an example of a `diag_current` section:

```text
diag_current
{
  ndump_fac = 20,
  reports = "j1", "j2", "j3" , 
}
```
