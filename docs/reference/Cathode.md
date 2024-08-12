---
layout: single
classes: wide
title: Cathode
permalink: /reference/cathode
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures cathode particle source settings. One of these
sections must exist for all `num_cathode` specified in the [particles](Particles.md)
section, and must be followed by a [species](Species.md) section, a (optional) [udist](Momentum_Distribution.md) section, a [spe_bound](Species_Boundary.md) section
and (optional) a [diag_species](Species_Diagnostics.md) section that refer to the species into which the
cathode particles will be injected. It accepts the following data:

- **dir**, integer, default = 1
- **wall**, integer, default = 1

<!-- -->

- **t_start**, real, default = 0.0
- **t_rise**, real, default = -1.0
- **t_fall**, real, default = -1.0
- **t_flat**, real, default = -1.0

<!-- -->

- **density**, real, default = 1.0
- **den_min**, real, default = 0.0
- **prof_type**, character(\*), default = "uniform"
- **mov_inj**, logical, default = .false.
- **deposit_current**, logical, default = .true.

<!-- -->

- **center**, real, default = -huge(1.0)
- **gauss_width**, real, default = -1.0
- **gauss_w0**, real, default = -1.0
- **channel_r0**, float, default = 0.0
- **channel_depth**, float, default = 0.0
- **channel_size**, float, default = 0.0
- **channel_wall**, float, default = 0.0
- **channel_bottom**, float, default = 1.0

**dir** - specifies the main direction along which to inject particles. It
must be in the range `[1,x_dim]`.

**wall** - specifies which wall to inject particles from. Valid values are
1 (lower wall) and 2 (upper wall).

**t_start**, **t_rise**, **t_fall**, **t_flat** - specify the temporal profile of the
particle bunch being injected. The `t_start` parameter specifies the time at which particle injection begins. The `t_rise` and `t_fall` parameters specify the particle
bunch rise and fall times, respectively. The shapes of the rise and fall are
approximately Gaussian and follow the usual 6th-order polynomial.
The `t_flat` parameter specifies the time the particle bunch holds its peak value. All
values are in normalized simulation units.

**density** - specifies a global multiplication factor for the transverse
density profile. Regardless of which profile type you choose, the final
density value will be `density * profile` value (see below). This parameter can be used
to set the peak density for the Gaussian profile and the density for the
uniform and step transverse density profile types.

**den_min** - specifies the minimum density for injecting particles.
Particles are only injected when the specified density is above this threshold.

**prof_type** - specifies the type of transverse profile to be used for
the particle bunch being injected. Valid values are:

- *"uniform"* - an uniform transverse profile
- *"gaussian"* - a Gaussian transverse profile
- *"step"* - a profile that is constant inside a given width and zero
  otherwise
- *"channel"* - a channel-like transverse profile (see parameters below)

**mov_inj** - specifies whether to use a static wall injection (default, `.false.`)
or a moving injector plane (`.true.`). Plane moves with the speed of light. Please
note that, to avoid unphysical fields, a neutral plasma must be injected
(zero currents).

**deposit_current** - specifies whether or not the electric current from the cathode particles should be deposited.

**center** - specifies the central position for the *"gaussian"* and *"step"*
transverse profiles. Note that when using either of these transverse
profiles, this parameter must be set.

**gauss_width** - specifies the full width to be used by the *"step"* profile,
and the maximum full width to be considered for the *"gaussian"* profile
(the value of the profile is set to zero if outside this maximum width).
Note that when using either of these transverse profiles, this parameter
must be set.

**gauss_w0** - specifies the waist (sigma) parameter of the *"gaussian"*
transverse profile. Note that when using the *"gaussian"* profile, this
parameter must be set.

**channel_r0**, **channel_depth**, **channel_size**,
**channel_wall**, **channel_bottom** - specify the parameters for a parabolic channel.
If `prof_type` is not set to *"channel"*, the
values specified are silently ignored. The density of the channel will
be `channel_bottom + channel_depth * (||x_perp -
channel_center||/channel_r0)^2` for `||x_perp - channel_center|| <
channel_size/2`, where `x_perp` is the position perpendicular to the
symmetry axis. For `||x_perp - channel_center|| >= channel_size/2`,
the density will either be constant (`channel_wall <= 0.0`) or fall off
linearly to 0.0 along a distance of `channel_wall`.

Note that any cathode particle source specification must be followed by
a [species](Species.md) section that specifies the parameters of the species into which the cathode
will inject particles (e.g., number of particles per cell, charge to mass ratio, etc.). Furthermore, the fluid velocity and thermal spread of the particles being injected are respectively set using the `ufl` and `uth`
parameters of an accompanying [udist](Momentum_Distribution.md) section. Finally, there must also a `spe_bound` section and an optional `diag_species` section, as with any other
species.

Here's an example of a cathode section followed by the definitions for
the associated species.

```text
cathode
{
  dir = 1,  ! direction to inject particles
  wall = 1, ! wall to inject particles from

  ! time profile information
  t_rise = 0.1,
  t_flat = 1000.0,
  t_fall = 3.0,

  ! transverse profile information
  density = 1.0,
  prof_type = "gaussian",
  center = 6.4,
  gauss_w0 = 3.2,
  gauss_width = 12.8,

  ! this sets the minimum density to inject particles from cathode
  den_min = 0.0000001d0,
}

!---------- new particles from cathode ----------
species
{
  name = "injected ions",

  num_par_max = 800000,
  rqm = -1836.0,
  num_par_x(1:2) = 2, 2,
}

!---------- initial proper velocities ----------
udist
{
  ! this sets the thermal momenta for the injected particles
  uth(1:3) = 7.2e-, 7.2e-, 7.2e-,

  ! this sets the fluid momenta for the injected   part.
  ufl(1:3)= 3.1e-, 0., 0.0d0,
}

!----------boundary conditions for this species----------
spe_bound
{
  type(1:2,1) = "open", "open",
  type(1:2,2) = "open", "open",
}

!----------diagnostic for this species----------
diag_species
{
  ndump_fac = 10,
  reports = "charge",
  ndump_fac_pha = 10,
  ndump_fac_raw = 10,
  ps_xmin(1:2) = 0.0, 0.0, ps_pmin(1:3) = -0.01, -0.01, -0.01,
  ps_xmax(1:2) = 12.8, 12.8, ps_pmax(1:3) = 0.01, 0.01, 0.01,
  ps_nx(1:2) = 128, 134, ps_np(1:3) = 256, 256, 256,
  if_ps_p_auto(1:3) = .true.,.true.,.true.,

  phasespaces = "p1x1", "p2x1", "p3x1", "p1x2", "p2x2", "p2p1",
}
```
