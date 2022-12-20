---
layout: single
classes: wide
title: Momentum Distribution
permalink: /reference/udist
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the momentum distribution for the particle
species and is optional. If not specified, the code will default to a
frozen (zero momentum) distribution. One of these sections should exist
for every species we intend to use. It accepts the following data:

- **uth_type**, character(\*), default = "thermal"

<!-- -->

- **uth**(3), real, default = 0.0
- **ufl**(3), real, default = 0.0
- **use_classical_uadd**, logical, default = .false.

<!-- -->

- **use_spatial_uth**, logical, default = .false.
- **use_spatial_ufl**, logical, default = .false.
- **spatial_uth**(3), character(\*), default = "0"
- **spatial_ufl**(3), character(\*), default = "0"

<!-- -->

- **math_func_uth**(3), character(\*), default = "0"
- **umin**(3), real, default = 0
- **umax**(3), real, default = 0
- **math_func_use_thermal**(3), logical, default = .false.

<!-- -->

- **relmax_T**, real, default = 0.0
- **relmax_umax**, real, default = -1.0

<!-- -->

- **n_accelerate**, integer, default = -1
- **n_q_incr**, integer, default = -1
- **n_accelerate_type**, character(\*), default = "linear"
- **n_q_incr_type**, character(\*), default = "linear"
- **use_particle_uacc**, logical, default = .false.

Particle momentum will be initialized as follows: i) the thermal
momentum is calculated in the particle rest frame and ii) this momentum
is boosted to the simulation frame using the fluid momentum.

**uth_type** - specifies the thermal distribution to use. Available
distributions are:

- *"none"* - Frozen (zero temperature) distribution
- *"thermal"* - Classical (Maxwell) exponential distribution
- *"random_dir"* - Random direction distribution, using module of uth
- *"waterbag"* - Waterbag distribution
- *"relmax"* - Relativistic (Juettner) exponential distribution
- *"waterbag_rel"* - Waterbag distribution with spherical distribution

**uth** - specifies the constant thermal spread in velocities for this
particle species in each of the directions. Momenta specified are proper
velocities, i.e., $\gamma \times v$ in units of c.

> For a *"thermal"* uth_type, the value given is the standard deviation of
a Maxwellian (i.e., Gaussian) and should be equal to $\sqrt{T/m c^2}$,
e.g., for a 1 keV isotropic electron plasma one would want to enter
`uth(1:3) = .04424, .04424, .04424,`.
>
> When specifying the thermal momenta for heavier species (i.e.,
$\left|rqm\right| > 1.0$), if we want them to have the same temperature
as an electron species, we need to set $u_{th}$ so that
$u_{th\,ion} = \sqrt{m_e / m_{ion}} u_{th\,e} \Leftrightarrow u_{th\,ion} = \sqrt{|1/rqm|} u_{th\,e}$.
For protons this would mean choosing
$u_{th} = u_{th\,e}/42.8504 = 0.0233370 \, u_{th\,e}$.
>
> When used in a species connected to a cathode particle source, it
specifies the thermal spread of the particles injected by the cathode.
This parameter is ignored if used in a species connected to a neutral or
neutral_mov_ions particle source.

**ufl** - specifies the constant fluid momentum for this particle species
in each of the directions. Momenta specified are proper velocities i.e.
$\gamma \times v$ in units of c. When used in a species connected to a
cathode particle source it specifies the fluid velocity of the particles
injected by the cathode. This parameter is ignored if used in a species
connected to a neutral or neutral_mov_ions particle source.

**use_classical_uadd** - setting this parameter to `.true.` forces a classical (linear) addition of
thermal velocity ($u_{th}$) and fluid velocity ($u_{fl}$). Relativistic
momentum addition is used by default.

**use_spatial_uth** - specifies that the thermal momentum is allowed to
have a spatial dependance. The thermal momenta are then defined using
the `spatial_uth` parameter. This is not compatible with the Juettner
distribution.

**use_spatial_ufl** - specifies that the fluid momentum is allowed to have
a spatial dependance. The fluid momenta are then defined using the
`spatial_ufl` parameter.

**spatial_uth** - specifies for each momentum component a
spatially dependent thermal momentum. To use this feature the user must
set the `use_spatial_uth` parameter. Momentum is specified using a
mathematical expression, which is a function of the spatial coordinates
`x1`, `x2` and `x3`. See the documentation on the analytical function parser
for details on the mathematical expression.

**spatial_ufl** - specifies for each momentum component a
spatially dependent fluid momentum. To use this feature the user must
set the `use_spatial_ufl` parameter. Momentum is specified using a
mathematical expression, which is a function of the spatial coordinates
`x1`, `x2` and `x3`. See the documentation on the analytical function parser
for details on the mathematical expression.

**math_func_uth** - specifies a function for each momentum component that is used to create a unique thermal distribution (e.g., something other than a Gaussian).  This mathematical expression is a function of `u`.  If a component is left unset or set to 0, a Gaussian will be used for that component.

**umin, umax** - specify the minimum and maximum allowed momenta for each momentum component when evaluating the `math_func_uth`.

**math_func_use_thermal** - specifies for each momentum component if the standard Gaussian should be used along with the associated thermal velocity.  This is akin to setting `math_func_uth = "0"` for a given direction.

**relmax_T** - specifies the temperature T for the Relativistic Maxwellian
distribution (relmax).

**relmax_umax** - specifies the maximum momentum to truncate the
Relativistic Maxwellian distribution. Please notice that the default
value (20 × relmax_T) was set only for convenience and may be too low
for high temperatures.

**n_accelerate** - specifies that the species is to be accelerated over
n_accelerate iterations along the x1 direction until it reaches the
fluid momentum specified in ufl(1). During acceleration particles are
only pushed using the longitudinal (x1) momentum. This is generally used
to initialize high energy particle beams with almost self-consistent
electro-magnetic fields. This parameter is ignored if used in a species
connected to a cathode, neutral or neutral_mov_ions particle source.

**n_q_incr** - specifies that the species is to be pushed over n_q_incr
iterations with increasing charge. During this time the particles are
free-streamed to initialize high energy particle beams with almost
self-consistent electro-magnetic fields. This parameter can be combined
with n_accelerate.

**n_accelerate_type**, **n_q_incr_type** - specifies kind of gradient to use
for the initialization of the species (see: `n_accelerate`, `n_q_incr`):

- *"linear"* - increase linearly
- *"regressive"* - increase regressive: first stronger, afterwards
  weaker (for example like: x/(x^2+1) or x/sqrt(x^2+1))

**use_particle_uacc** - setting this parameter `.true.` in combination with a particle source from file will accelerate particles to their $p_1$ as given in the RAW data file.  See extra documentation for initializing from a RAW file.

Here's an example of an udist section using a Maxwellian distribution,
where the temperature varies spatially, and with a fixed fluid momentum
of 10.0 in the x1 direction:

```text
udist  {  
  ! Non-uniform thermal  
  use_spatial_uth = .true.,  
  spatial_uth(1) = 'x1',  
  spatial_uth(2) = 'x2',  
  spatial_uth(3) = 'x1*x2',  

  ! Uniform fluid  
  ufl(1:3) =  10.0, 0.0, 0.0, 
}
```
