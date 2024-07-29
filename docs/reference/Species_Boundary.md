---
layout: single
classes: wide
title: Species Boundary
permalink: /reference/spe_bound
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures particle species boundary confition setting and
must be present in the input file. One of these sections must exist for
every species we intend to use. It accepts the following data:

- **type**(2,x_dim), character(\*), default = "---"
- **uth_bnd**(3,2,x_dim), float, default = 0.0
- **ufl_bnd**(3,2,x_dim), float, default = 0.0
- **thermal_type**(2,x_dim), character(\*), default = "thermal"

**type** - specifies the boundary conditions to use for the particle
species. The `type(k,i)`, with `k = {1,2}` and `i = {1, ..., x_dim}` specifies the
boundary condition in the `i` direction for the lower (`k=1`) or upper (`k=2`)
boundary. The following are valid values for type:

- "axial" - Axial boundary condition for 2D cylindrical coordinates.
- "absorbing" or "open" - Absorbing boundary conditions.
- "specular" or "reflecting" - Specular reflection boundary conditions.
- "thermal" - Thermal bath boundary conditions. See `uth_bnd` and
  `ufl_bnd` for details.

If periodic boundaries or a moving window has been specified for a given
direction, then `type` is ignored for that direction.

**uth_bnd**, **ufl_bnd** - specify the parameters for the thermal bath
boundary conditions. When using a thermal bath boundary condition, any
particle crossing the specified boundary is re-injected into the
simulation at its previous position using the specified thermal spread
(`uth_bnd`) and fluid velocity (`ufl_bnd`). Like the velocity parameters
in the species section, these are specified as proper velocities, i.e.,
$\gamma \times v$ in units of c, and the same comments apply. The `uth_bnd(j,k,i)`
and `ufl_bnd(j,k,i)`, with `j = {1..3}`, `k = {1,2}` and `i = {1, ..., x_dim}`
specify the `j` component of the thermal and fluid velocities,
respectively, in the `i` direction for the lower (`k=1`) or upper (`k=2`)
boundary. if `type(k,i)` is not set to "thermal", the values specified are
silently ignored.

**thermal_type** - specifies the type of thermal distribution to be used
for thermal boundary conditions. The `thermal_type(k,i)`, with `k = {1,2}` and `i
= {1, ..., x_dim}` specify the thermal distribution type in the `i`
direction for the lower (`k=1)` or upper (`k=2`) boundary. Possible values
are:

- _"thermal"_ - Use a half-Maxwellian distribution for the velocity
  component perpendicular to the wall, and use the standard classical
  (Maxwellian) exponential distribution for all other components.
- _"waterbag"_ - Waterbag distribution.
- _"random dir"_ - Random direction distribution, using absolute value
  of `uth`.

In all cases, the resulting velocity component perpendicular to the wall
will point inwards.

Here's an example of an spe_bound section for a 2D run that sets the the
boundary conditions in the x1 direction to absorbing boundaries and in
the x2 direction to thermal bath boundaries. The thermal velocities for
the boundaries in the x2 direction are set to (0.1,0.1,0.1) on the lower
boundary and to (0.001,0.001,0.001) on the upper boundary.

```text
  spe_bound
  {
    type(1:2,1) =    "open",    "open",
    type(1:2,2) = "thermal", "thermal",
    uth_bnd(1:3,1,2) = 0.1,0.1,0.1,
    uth_bnd(1:3,2,2) = 0.001,0.001,0.001,
  }
```

Note: When using thermal boundary conditions, it is necessary to modify the
settings for the field boundary conditions from the default. In 2D cartesian
geometry, the option **vpml_diffuse** enables damping of noise from absorbing
particles that cross the simulation boundary.

```text
    emf_bound
    {
      type(1:2,1) = "vpml", "vpml",         ! Equivalent to "open", but more explicit
      type(1:2,2) = "vpml", "vpml",
      vpml_bnd_size = 30,                   ! Found to provide good results in a 2D test, larger than the default of 8
      vpml_diffuse(1:2,1) = .true., .true., ! Allows for damping of electrostatic fields in vpml
      vpml_diffuse(1:2,2) = .true., .true.,
    }
```
