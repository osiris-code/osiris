# Neutrals with Moving Ions

This section configures neutral with moving ions settings and must be
present in the input file if any neutrals with moving ions were
specified with the parameter num_neutral_mov_ions in the particles
section. This module injects electrons from a tunnel ionizable fixed gas
background with moving ions (see also the neutral section). One of these
sections must exist for every neutral with moving ions we intend to use,
and must be followed by two sets of species, spe_bound and (optional)
diag_species sections that refer to the species the ionization electrons
and moving ions will be injected into. It accepts the following data:

- **name**, string, default = "Neutral"
- **neutral_gas**, string, default = "H"
- **ion_param**(3,ion_max), float, default = 0.0
- **den_min**, float, default = 0.0
- **e_min**, float, default = 1.0e-6
- **multi_max**, integer, default = maximum ionization level available
- **multi_min**, integer, default = 0
- **inject_line**, logical, default = .true.

The associated species must define the number of particles per cell that
will be used. The thermal and fluid velocities specified will be
ignored. This module can optionally be followed by a profile section (or
one num_x and one profile section), placed before the species
information, that will define the density of the neutral background gas
density.

**neutral_name** specifies the name for the neutral.

**neutral_gas** specifies the neutral gas to use in the simulation.
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
- *"custom"* - Use a custom ionization model defined by the *ion_param*
  parameter.

**ion_param** Allows the user to specify a custom neutral background
gas. The ionization rate for each ionization level is given in $s^{-1}$
by \[Bruhweiler D et al., Physics of Plasmas 10, 2022 (2003), eqs. (1)
and (2)\] :

$W[s^{-1}] = A \times E^{-C} \times e^{-B/E}$

Where $E$ is the absolute value of the electric field in GV/m, and the
$A$, $B$ and $C$ parameters are defined for ionization level i as *A* =
ion_param(1,i), *B* = ion_param(2,i), *C* = ion_param(3,i). See the
above mentioned paper for details on calculating these values for an
arbitrary medium.

**den_min** specifies the minimal density of neutrals to be considered
for injection.

**e_min** specifies the minimum field value, in GV/m to be considered
for ionization. If the field is below this value it will not ionize the
gas.

**multi_max** specifies the maximum number of ionization levels. It can
be used to truncate ionization at a given level.

**multi_min** specifies the minimum ionization level after which
particles are injected (useful to assign a different species to a set of
levels).

**inject_line** specifies whether to inject particles uniformly along a
line in the x1 direction (.true.) or to place them randomly in the cell
(.false.)

Here's an example of a neutral_mov_ions section using a Carbon
background and allowing ionization only up to the 3rd level:

```text
neutral_mov_ions { 
  neutral_gas = "C",
  multi_max = 3,
}
```
