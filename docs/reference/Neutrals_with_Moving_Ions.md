---
layout: single
classes: wide
title: Neutrals with Moving Ions
permalink: /reference/neutral_mov_ions
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the settings of neutrals with moving ions and must be
present in the input file if any neutrals with moving ions were
specified with the parameter `num_neutral_mov_ions` in the [particles](Particles.md)
section. This module injects electrons from a tunnel-ionizable fixed gas
background with moving ions (see also the [neutral](Neutrals.md) section). One of these
sections must exist for every neutral with moving ions we intend to use,
and must be followed by two sets of [species](Species.md), [spe_bound](Species_Boundary.md) and (optional)
[diag_species](Species_Diagnostics.md) sections that refer to the species into which the ionization electrons
and moving ions will be injected, respectively. It accepts the following data:

- **name**, string, default = "Neutral"
- **neutral_gas**, string, default = "H"
- **ion_param**(3,ion_max), float, default = 0.0
- **den_min**, float, default = 0.0
- **e_min**, float, default = 1.0e-6
- **multi_max**, integer, default = maximum ionization level available
- **multi_min**, integer, default = 0
- **inject_line**, logical, default = .true.
- **if_tunnel**, logical, default = .true.
- **if_impact**, logical, default = .false.

The associated species must define the number of particles per cell that
will be used. The thermal and fluid velocities specified will be
ignored. This module can optionally be followed by a [profile](Profile.md) section,, placed before the [species](Species.md)
section, that will define the density of the neutral background gas
density.

**name** - specifies the name for the neutral.

**neutral_gas** - specifies the neutral gas to use in the simulation.
Available gases are:

- *"H"* - Hydrogen
- *"Li"* - Lithium
- *"Na"* - Sodium
- *"K"* - Potassium
- *"Rb"* - Rubidium
- *"Cs"* - Cesium
- *"He"* - Helium
- *"Ar"* - Argon
- *"N"* - Nitrogen
- *"O"* - Oxygen
- *"C"* - Carbon
- *"Ne"* - Neon
- *"Xe"* - Xenon
- *"custom"* - Use a custom ionization model defined by the `ion_param`
  parameter.

**ion_param** - Allows the user to specify a custom neutral background
gas. The ionization rate for each ionization level is given in $s^{-1}$
by \[Bruhweiler D et al., Physics of Plasmas 10, 2022 (2003), eqs. (1)
and (2)\]:

$W[s^{-1}] = A \times E^{-C} \times e^{-B/E}$

Where $E$ is the absolute value of the electric field in GV/m, and the
$A$, $B$ and $C$ parameters are defined for ionization level i as `A =
ion_param(1,i)`, `B = ion_param(2,i)` and `C = ion_param(3,i)`. See the
above mentioned paper for details on calculating these values for an
arbitrary medium.

**den_min** - specifies the minimal density of neutrals to be considered
for injection.

**e_min** - specifies the minimum field value, in GV/m to be considered
for ionization. If the field is below this value, it will not ionize the
gas.

**multi_max** - specifies the maximum number of ionization levels. It can
be used to truncate ionization at a given level.

**multi_min** - specifies the minimum ionization level after which
particles are injected (useful to assign a different species to a set of
levels).  Keyword should be set to the desired minimum ionization level
minus one.

**inject_line** - specifies whether to inject particles uniformly along a
line in the `x1` direction (`.true.`) or to place them randomly in the cell
(`.false.`).

**if_tunnel** - specifies whether or not tunnel ionization is turned on.

**if_impact** - specifies whether or not impact ionization is turned on.

Here's an example of a neutral_mov_ions section using a Carbon
background and allowing ionization only up to the 3rd level:

```text
neutral_mov_ions
{
  neutral_gas = "C",
  multi_max = 3,
}
```
