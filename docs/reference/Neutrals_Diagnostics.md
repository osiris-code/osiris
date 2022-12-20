---
layout: single
classes: wide
title: Neutrals Diagnostics
permalink: /reference/diag_neutral
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the neutral diagnostic settings and is optional.
If not present, the code will not do any neutral diagnostics. It accepts the following data:

- **reports**(:), character(\*), default = "-"
- **ndump_fac**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **n_tavg**, integer, default = -1
- **n_ave**(x_dim), integer, default = -1
- **prec**, integer, default = 4

**reports** - specifies the grid quantities to report, including
spatial/time averaging, lineouts, etc., as described in the [grid diagnostics section](../other/grid_diagnostics). The
available quantities are:

- "ion_charge" - Background ion charge density (particle density times
  ionization level)
- "neut_den" - Density of the initial background neutral gas

**ndump_fac** - controls the frequency of full grid diagnostics. This
value is multiplied by the `ndump` value specified in the [time_step](Time_Step.md)
section to determine the number of iterations between each diagnostic
dump. If set to 0, the writing of this diagnostic information is
disabled.

**ndump_fac_ave** - controls the frequency of spatial average / envelope
grid diagnostics. This value is multiplied by the `ndump` value
specified in the [time_step](Time_Step.md) section to determine the number of
iterations between each diagnostic dump. If set to 0, the writing of this
diagnostic information is disabled.

**ndump_fac_lineout** - controls the frequency of lineout / slice
diagnostics. This value is multiplied by the `ndump` value specified in
the [time_step](Time_Step.md) section to determine the number of iterations between
each diagnostic dump. If set to 0, the writing of this diagnostic
information is disabled.

**n_tavg** - specifies the number of time steps to be used when
calculating the time averaged diagnostics. The frequency of these
diagnostics is controlled by the `ndump_fac` parameter described above.

**n_ave** - number of gridpoints on each direction to average over for
spatially averaged dumps. The frequency of these diagnostics is
controlled by the `ndump_fac_ave` parameter described above.

**prec** - controls the numerical precision used for grid diagnostics. The
default is to save data in single precision (`prec = 4`). If the user
wants data to be saved in double precision, this parameter must be set to
8. This option is ignored if OSIRIS is compiled in single precision.

Here's an example of a diag_neutral section that will write diagnostics
information every `20*ndump` iterations for ion charge density.

```text
diag_neutral
{
  ndump_fac = 20,
  reports = "ion_charge",
}
```
