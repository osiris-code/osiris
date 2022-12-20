---
layout: single
classes: wide
title: Species Diagnostics
permalink: /reference/diag_species
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the particle species diagnostic settings and is
optional. If not present, the code will not do any species diagnostics.
It accepts the following data:

- **reports**(:), character(\*), default = "-"
- **rep_cell_avg**(:), character(\*), default = "-"
- **rep_udist**(:), character(\*), default = "-"
- **ndump_fac**, integer, default = 0
- **ndump_fac_ene**, integer, default = 0
- **ndump_fac_temp**, integer, default = 0
- **ndump_fac_heatflux**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **n_tavg** , integer, default = 0
- **n_ave**(x_dim), integer, default =0
- **prec**, integer, default = 4

<!-- -->

- **ndump_fac_raw**, integer, default = 0
- **raw_gamma_limit**, real, default = 0.0
- **raw_fraction**, real, default = 1.0
- **raw_math_expr**, character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **raw_if_pos_ref_box**(x_dim), logical, default = .false.

<!-- -->

- **ndump_fac_tracks**, integer, default = 0
- **n_start_tracks**, integer, default = -1
- **niter_tracks**, integer, default = 1
- **file_tags**, character(\*), default = ""
- **ifdmp_tracks_efl**(3), bool, default = .false.
- **ifdmp_tracks_bfl**(3), bool, default = .false.
- **ifdmp_tracks_psi**, bool, default = .false.

<!-- -->

- **phasespaces**(:), character(\*), default = "-"
- **pha_ene_bin**(:), character(\*), default = "-"
- **pha_cell_avg**(:), character(\*), default = "-"
- **pha_time_avg**(:), character(\*), default = "-"
- **ndump_fac_pha**, integer, default = 0
- **ndump_fac_pha_tavg**, integer, default = 0

<!-- -->

- **ps_xmin**(x_dim), real, default = 0.0
- **ps_xmax**(x_dim), real, default = 0.0
- **ps_nx**(x_dim), integer, default = 64
- **ps_nx_3D**(x_dim), integer, default = ps_nx

<!-- -->

- **ps_pmin**(3), real, default = 0.0
- **ps_pmax**(3), real, default = 0.0
- **ps_np**(3), integer, default = 64
- **ps_np_3D**(3), integer, default = ps_np
- **if_ps_p_auto**(3), bool, default = .false.

<!-- -->

- **ps_lmin**(3), real, default = 0.0
- **ps_lmax**(3), real, default = 0.0
- **ps_nl**(3), integer, default = 64
- **ps_nl_3D**(3), integer, default = ps_np
- **if_ps_l_auto**(3), bool, default = .false.

<!-- -->

- **ps_gammamin**, real, default = 1.0
- **ps_gammamax**, real, default = 0.0
- **ps_ngamma**, integer, default = 64
- **if_ps_gamma_auto**, bool, default = .false.

<!-- -->

- **ps_kemin**, real, default = 0.000001
- **ps_kemax**, real, default = 0.0
- **ps_nke**, integer, default = 64
- **if_ps_ke_auto**, bool, default = .false.

<!-- -->

- **n_ene_bins**, integer, default = 0
- **ene_bins**(:), real, default = 0

**reports** - specifies the spatially resolved quantities to report,
including spatial/time averaging, lineouts, etc., as described in the
[grid diagnostics section](../other/grid_diagnostics). The available
quantities are:

- *"charge"* - Species charge
- *"m"* - Species mass ( charge\*rqm )
- *"ene"* - Local kinetic energy
- *"q1"*, *"q2"*, *"q3"* - Heat flux component
- *"j1"*, *"j2"*, *"j3"* - Electric current component for this species. Please
  note that this current is not exactly the same (but very close) as the
  current used in the OSIRIS algorithm. This current is calculated
  depositing q\*v on a grid, and the OSIRIS current is calculated using a
  charge-conserving scheme. The differences between the two are, however,
  minimal.

**rep_cell_avg** - specifies the cell averaged quantities to report,
including spatial/time averaging, lineouts, etc., as described in the
[grid diagnostics section](../other/grid_diagnostics). The available
quantities are the same as for the `reports` item. After depositing the
selected quantity, it will be divided by the absolute cell density.
These diagnostics are controlled by the same parameters as the ones
specified by the `report` item.

**rep_udist** - specifies grid resolved fluid and thermal momenta
diagnostics, including spatial/time averaging, lineouts, etc., as
described in the [grid diagnostics section](../other/grid_diagnostics). The available
quantities are:

- *"ufl1"*, *"ufl2"*, *"ufl3"* - Average fluid momentum component.
- *"uth1"*, *"uth2"*, *"uth3"* - Momentum distribution width along the
  specified direction.

**ndump_fac** - controls the frequency of full grid diagnostics. This
value is multiplied by the `ndump` value specified in the [time_step](Time_Step.md)
section to determine the number of iterations between each diagnostic
dump. If set to 0, the writing of this diagnostic information is
disabled.

**ndump_fac_ene** - specifies the frequency at which to write total
particles species energy diagnostics. This value is
multiplied by the `ndump` value specified in the [time_step](Time_Step.md) section to
determine the number of iterations between each diagnostic dump. If set
to 0, the writing of this diagnostic information is disabled.

**ndump_fac_temp** - specifies the frequency at which to write total
particles species temperature diagnostics. This value is multiplied by
the `ndump` value specified in the [time_step](Time_Step.md) section to determine the
number of iterations between each diagnostic dump. If set to 0, the
writing of this diagnostic information is disabled.

**ndump_fac_heatflux** - specifies the frequency at which to write total
particles species heat flux diagnostics. This value is multiplied by
the `ndump` value specified in the [time_step](Time_Step.md) section to determine the
number of iterations between each diagnostic dump. If set to 0, the
writing of this diagnostic information is disabled.

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

### Raw Diagnostics

**ndump_fac_raw** - specifies the frequency at which to write particle
species raw diagnostics. This diagnostic dumps all particle information
(position, momenta and charge) to file. This value is multiplied by the
`ndump` value specified in the [time_step](Time_Step.md) section to determine the number
of iterations between each diagnostic dump. If set to 0, the writing of
this diagnostic information is disabled. See also `raw_gamma_limit` and
`raw_fraction`.

**raw_gamma_limit** ( \>= 1.0 ) - minimal relativistic gamma for raw
diagnostic. Only particle data from particles with `gamma >= gamma_limit`
will be saved. Note that this parameter has nothing to do with the gamma
distribution diagnostic.

**raw_fraction** ( \[0.0, 1.0\] ) - fraction of particles to dump for raw
diagnostic. Particles are selected by randomly testing for `random <=
particle_fraction` so that approximately only particle data from
`particle_fraction` of the total particles are saved.

**raw_math_expr** - specifies a mathematical function for particle
selection in the raw diagnostic. This expression can be a function of
any of the following:

- **x{1\|2\|3}** - represent the physical coordinates of the particle.
- **p{1\|2\|3}** - represent the generalized momenta of the particle.
- **g** - represents the particle Lorentz gamma factor.
- **t** - represents the current simulation time.

See the documentation on the analytical function parser for details on
the mathematical expression.

Note that `raw_gamma_limit`, `raw_fraction` and `raw_math_expr` may all be specified and actively used to determine which particles are selected for the raw diagnostic

**raw_if_pos_ref_box** - specify in each direction whether the particle position in the raw diagnostic will be saved as its physical position (`.false.`) or in reference to the box edge (`.true.`). This is most often set to true in the moving window direction to prevent roundoff errors.

### Particle tracking

**ndump_fac_tracks** - specifies the frequency at which to write particle
tracking information to file. This value is multiplied by the `ndump`
value specified in the [time_step](Time_Step.md) section to determine the number of
iterations between each diagnostic dump. All the values for particles
whose tags are in the `file_tags` file will be saved to memory every
`niter_tracks` time steps. This information will be flushed to disk at
a frequency defined by this parameter. Note that for this diagnostic, you
also need to set `add_tags` in the parent [species](Species.md) section.

**n_start_tracks** - specifies the iteration after which tracking information will be collected.

**niter_tracks** - specifies the number of iterations between saving a
track point. For example, setting this parameter to 10 would mean that
track points would be added at every 10 time steps.

**file_tags** - specifies the name of the file that holds the list of tags
of the particles we want to follow.

**ifdmp_tracks_efl** - specifies whether to save the electric field in the
tracks file for each direction

**ifdmp_tracks_bfl** - specifies whether to save the magnetic field in the
tracks file for each direction

**ifdmp_tracks_psi** - specifies whether to save the psi quantity in the
tracks file for each direction

### Phasespaces

**phasespaces** - specifies the phasespace diangnostics to perform. This
parameter is an array of strings, where each string specifies the
parameters for each of the diagnostics to perform. The strings specify
which quantities to use in the phasespace simply by using the quantity
name, e.g., `"x1"` will do a 1D phasespace with the `x1` position, `"p1x1"` will
do a standard 2D momentum/position (`p1`/`x1`) phasespace and `"p3x2x1"` will
do a 3D phasespace with `p3`, `x2` and `x1`. Valid quantities for use as
phasespace axis are as follows:

- **x{1\|2\|3}** – particle position along directions 1, 2 and 3
  respectively. Specifying `x3` in a 2D run or specifying `x3` or `x2` in a 1D
  run results in a run-time error.
- **p{1\|2\|3}** – particle linear momentum along directions 1, 2 and 3
  respectively.
- **l{1\|2\|3}** – particle angular momentum along directions 1, 2 and 3
  respectively, calculated referenced to the center of the simulation
  box. If using a moving window, the angular momentum is calculated
  referencing the center of the current simulation window.
- **g**, **gl** - lorentz factor of the particle and logarithm of the lorentz
  factor of the particle, respectively.

By default, the charge of the particle is used in the phasespace, but
other quantities can be used by appending “\_QUANT” to the phase-space
string, where QUANT stands for the quantity being deposited. Quantities
that can be deposited are the following:

- **charge** (default) – deposits the charge of the particle, ex: `“x1”` or
  `“x1_charge”`.
- **m** – deposits the mass of the particle, ex: `“x2x1_m”`.
- **ene** – deposits the kinetic energy of the particle, ex: `“x1_ene”`.
- **\|charge\|** - deposits the absolute value of the charge of the
  particle, ex: `“x3x2x1_|charge|”`.
- **q{1\|2\|3}** – deposits the heat flux of the particle in the given
  direction, ex: `“x2x1_q2”` would deposit the heat flux along `x2` in an
  `x2x1` phasespace.
- **j{1\|2\|3}** – deposits the electrical current of the particle in the
  given direction, ex: `“x3x1_j1”` would deposit the current along `x1` in an
  `x3x1` phasespace.

The maximum number of phasespaces that a user can specify is defined by
the constant `p_max_phasespaces`, set at compile time in os-spec-diagnostics.f03. The
default is 64.

**pha_ene_bin** - specifies the energy binned phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
is separated into a set of phasespaces that only take into account the
particles with energies falling into the corresponding energy bin. The
number of energy bins and the ranges for each bin are specified through
the `n_ene_bins` and `ene_bins` parameters.

**pha_cell_avg** - specifies the cell averaged phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
stores the average value of the quantity being deposited for every
particle inside each cell (rather than the total sum of this quantity).
If no particles are present inside a given cell, the value in that cell
is set to zero.

**pha_time_avg** - specifies the time averaged phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
stores the average value of the phasespace over a number of time steps
specified by the `n_tavg` parameter. Please note that 1) it is not
possible to make time averaged phasespaces if any of the axis used are
being autoranged, as the ranges may change while data is being averaged,
and 2) for each time averaged phasespace that the user requests, the
code will need to hold its grid in memory, which means that using these
phasespaces increases the memory requirements for OSIRIS significantly.

**ndump_fac_pha** - specifies the frequency at which to write phasespace
particle species diagnostics. This value is multiplied by the `ndump`
value specified in the [time_step](Time_Step.md) section to determine the number of
iterations between each diagnostic dump. If set to 0, the writing of this
diagnostic information is disabled. The individual phasespace dumps are
turned on and off through the phasespaces, pha_ene_bin and pha_cell_avg
parameters. See these parameters for details.

**ndump_fac_pha_tavg** - specifies the frequency at which to write time
averaged phasespace particle species diagnostics. This value is
multiplied by the `ndump` value specified in the [time_step](Time_Step.md) section to
determine the number of iterations between each diagnostic dump. If set
to 0, the writing of this diagnostic information is disabled. The
individual phasespace dumps are turned on and off through the
`pha_time_avg` parameter. See this parameter for details.

**ps_xmin**, **ps_xmax** - specify the lower and upper limit, for every
direction, to be considered for phasespace diagnostics using the `x1`, `x2`
or `x3` (in 3D) coordinate.

**ps_nx** - specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the `x1`, `x2` or `x3`
(in 3D) coordinate.

**ps_nx_3D** - specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the `x1`, `x2` or `x3` (in 3D).
If not specified, it defaults to the same value as `ps_nx`. This parameter
allows for these 3D phasespaces to be done at lower resolutions while
maintaining a high resolution for 1D and 2D phasespaces.

**ps_pmin**, **ps_pmax** - specify the lower and upper limit, for every
direction, to be considered for 1D and 2D phasespace diagnostics using
the `p1`, `p2` or `p3` quantities (linear momentum).

**ps_np** - specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the `p1`, `p2` or `p3`
quantities.

**ps_np_3D** - specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the the `p1`, `p2` or `p3`
quantities. If not specified defaults to the same value as `ps_np`. This
parameter allows for these 3D phasespaces to be done at lower
resolutions while maintaining a high resolution for 1D and 2D
phasespaces.

**if_ps_p_auto** - specifies whether the code should try to autoscale the
momenta phasespace limits set by `ps_pmin` and `ps_pmax` so that all
particles in the species are present in the phasespace diagnostic. The
code will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by `ps_pmin` and `ps_pmax` will always be present
in the phasespace (the code will never shrink below these parameters).
Note that you can select autoscaling in each momentum direction
independently, i.e., each component of `if_ps_p_auto` corresponds to
autoscaling the momentum in that direction.

**ps_lmin**, **ps_lmax** - specify the lower and upper limit, for every
direction, to be considered for 1D and 2D phasespace diagnostics using
the `l1`, `l2` or `l3` quantities (angular momentum).

**ps_nl** - specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the `l1`, `l2` or `l3`
quantities.

**ps_nl_3D** - specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the the `l1`, `l2` or `l3`
quantities. If not specified defaults to the same value as `ps_nl`. This
parameter allows for these 3D phasespaces to be done at lower
resolutions while maintaining a high resolution for 1D and 2D
phasespaces.

**if_ps_l_auto** - specifies whether the code should try to autoscale the
momenta phasespace limits set by `ps_lmin` and `ps_lmax` so that all
particles in the species are present in the phasespace diagnostic. The
code will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by `ps_lmin` and `ps_lmax` will always be present
in the phasespace (the code will never shrink below these parameters).
Note that you can select autoscaling in each momentum direction
independently, i.e., each component of `if_ps_l_auto` corresponds to
autoscaling the angular momentum in that direction.

**ps_gammamin**,**ps_gammamax** - specify the lower and upper limit to be
considered for the gamma distribution diagnostics.

**ps_ngamma** - specifies the number of points to be
considered for 1D, 2D or 3D phasespace diagnostics using gamma.

**if_ps_gamma_auto** - specifies whether the code should try to autoscale
the gamma distribution limits set by `ps_gammamin` and `ps_gammamax` so that
all particles in the species are present in this diagnostic. The code
will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by `ps_gammamin` and `ps_gammamax` will always be
present in the diagnostic (the code will never shrink below these
parameters).

**ps_kemin**,**ps_kemax** - specify the lower and upper limit to be
considered for the kinetic energy distribution diagnostics.

**ps_nke** - specifies the number of points to be
considered for 1D, 2D or 3D phasespace diagnostics using kinetic energy.

**if_ps_ke_auto** - specifies whether the code should try to autoscale
the kinetic energy distribution limits set by `ps_kemin` and `ps_kemax` so that
all particles in the species are present in this diagnostic. The code
will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by `ps_kemin` and `ps_kemax` will always be
present in the diagnostic (the code will never shrink below these
parameters).

**n_ene_bins** - specifies the number of values that the `ene_bins`
parameter holds. The actual number of energy bins used for the energy
binned phasespace diagnostics defined by the `pha_ene_bin` parameter are
actually `n_ene_bins + 1`. See `ene_bins` and `pha_ene_bins` for details.

**ene_bins** - specifies the energies defining the energy bins to be used
by the energy binned phasespace diagnostics defined by the `pha_ene_bin`
parameter. The energy for each particle is calculated as `gamma - 1`
(particle kinetic energy normalized to rest mass), where gamma is the
relativistic Lorentz factor. The first bin will hold the particles with
`energies <= ene_bins(1)`, the last bin will hold particles with
`energies > ene_bins(n_ene_bins)` and the i-th bin will hold particles
with `ene_bin(i-1) < energy <= ene_bin(i)`. The results for all bins are
saved in the same file.

Here's an example of a `diag_species` section:

```text
diag_species
{
  ndump_fac = 1,
  reports = "charge",

  ndump_fac_pha = 1,
  ndump_fac_raw = 1,
  ps_xmin(1:2) =  0.0, 0.0,
  ps_xmax(1:2) = 15.0, 3.5,
  ps_nx(1:2)   = 4096, 512,

  ps_pmin(1:3) = -2.0, -2.0, -2.0,
  ps_pmax(1:3) = 40.0,  2.0,  2.0,
  ps_np(1:3)   = 1024,  256,  256,
  if_ps_p_auto(1:3) = .true., .true., .true.,

  ps_gammamin = 1.0,
  ps_gammamax = 50.0,
  ps_ngamma = 8192,
  if_ps_gamma_auto = .true.,

  phasespaces = "x1", "gl_m", "p1x1", "x2x1_ene", "p2x1x2",

  raw_gamma_limit = 10.0,
  raw_fraction = 1.0,
}
```
