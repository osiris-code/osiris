---
layout: single
classes: wide
title: Particles
permalink: /reference/particles
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the particle sources to be used in the
simulation. These consist mainly of two types: particles injected at the
begining of the simulation obeying some spatial profile and temperature
distribution ([species](Species.md)), and
sources that inject particles into the simulation during the run, either
from some external source
([cathode](Cathode.md)) or from ionization
([neutral](Neutrals.md),
[neutral_mov_ions](Neutrals_with_Moving_Ions.md)).

- **num_species**, integer, default = 0
- **num_cathode**, integer, default = 0
- **num_neutral**, integer, default = 0
- **num_neutral_mov_ions**, integer, default = 0

<!-- -->

- **interpolation**, character(\*), default = "quadratic",
- **low_jay_roundoff**, logical, default = .false.
- **grid_center**, logical, default = .false.

<!-- -->

- **reports**(:), character(\*), default = "-"
- **ndump_fac**, integer, default = 0
- **ndump_fac_ene**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **n_tavg**, integer, default = -1
- **n_ave**(x_dim), integer, default = -1
- **prec**, integer, default = 4

**num_species** - specifies the number of particle species injected at the
begining of the simulation to use. See also the
[species](Species.md) section for details.

**num_cathode** - specifies the number of cathodes to use in the
simulation. See also the [cathode](Cathode.md)
section for details.

**num_neutral** - specifies the number of neutral (tunnel ionizable gas)
to use in the simulation. See also the
[neutral](Neutrals.md) section for details.

**num_neutral_mov_ions** - specifies the number of neutral_mov_ions
(tunnel ionizable gas with moving ions) to use in the simulation. See
also the
[neutral_mov_ions](Neutrals_with_Moving_Ions.md)
section for details.

**interpolation** - specifies the interpolation level to be used for all
particles. Possible values are:

- *"linear"* - Use linear interpolation
- *"quadratic"* - Use quadratic interpolation
- *"cubic"* - Use cubic interpolation
- *"quartic"* - Use quartic interpolation

**low_jay_roundoff** - controls the way current from multiple species is added
to obtain the total current. The default behavior is to sequentially add
current from multiple species directly to the global current grid.
However, for cases with zero net current (at initialization), the
differences due to roundoff will cause a numeric noise of the order
$10^{-12} - 10^{-15}$ that can lead to instabilities.
Setting this parameter to `.true.` will deposit each species on a clean
grid, then add the grid results from all species. This technique
effectively eliminates the numeric noise (roundoff errors are
independent of the sign), but has the caveat of requiring more memory for
the extra current grids.

**grid_center** - specifies if the field values are centered at the corner of the cell (full momentum conserving).

**reports** - specifies the grid quantities to report, including
spatial/time averaging, lineouts, etc., as described in the [grid diagnostics section](../other/grid_diagnostics). The
available quantities are:

- "charge" - Total simulation charge
- "charge_htc" - Total simulation charge, time centered at n-½ (i.e.,
  same as electric current)
- "dcharge_dt" - Spatially resolved time derivative of total simulation
  charge

**ndump_fac** - controls the frequency of full grid diagnostics. This
value is multiplied by the `ndump` value specified in the `time_step`
section to determine the number of iterations between each diagnostic
dump. If set to 0, the writing of this diagnostic information is
disabled.

**ndump_fac_ene** - controls the frequency of energy diagnostics ofr all species. This
value is multiplied by the `ndump` value specified in the `time_step`
section to determine the number of iterations between each diagnostic
dump. If set to 0, the writing of this diagnostic information is
disabled.

**ndump_fac_ave** - controls the frequency of spatial average / envelope
grid diagnostics. This value is multiplied by the `ndump` value
specified in the `time_step` section to determine the number of
iterations between each diagnostic dump. If set to 0, the writing of this
diagnostic information is disabled.

**ndump_fac_lineout** - controls the frequency of lineout / slice
diagnostics. This value is multiplied by the `ndump` value specified in
the `time_step` section to determine the number of iterations between
each diagnostic dump. If set to 0, the writing of this diagnostic
information is disabled.

**n_tavg** - specifies the number of time steps to be used when
calculating the time-averaged diagnostics. The frequency of these
diagnostics is controlled by the `ndump_fac` parameter described above.

**n_ave** - number of gridpoints in each direction to average over for
spatially averaged dumps. The frequency of these diagnostics is
controlled by the `ndump_fac_ave` parameter described above.

**prec** - controls the numerical precision used for grid diagnostics. The
default is to save data in single precision (`prec = 4`). If the user
wants data to be saved in double precision, this parameter must be set to 8. This option is ignored if OSIRIS is compiled in single precision.

This section must be followed by `num_species`
[species](Species.md) sections, `num_cathode`
[cathode](Cathode.md) sections, `num_neutral`
[neutral](Neutrals.md) sections, and `num_neutral_mov_ions`
[neutral_mov_ions](Neutrals_with_Moving_Ions.md)
sections, in this order, describing each individual particle source.

Here's an example of a particles section setting the number of normal
species to 2, and the number of cathodes to 1, and setting the
interpolation to linear interpolation.

```text
particles
{
  num_species = 2,
  num_cathode = 1,

  interpolation = "linear",
}
```
