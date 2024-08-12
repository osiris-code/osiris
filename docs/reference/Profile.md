---
layout: single
classes: wide
title: Profile
permalink: /reference/profile
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the density profile for the particle species and
is generally required in the input file. If parameters are not specified, the
code will default to a uniform density profile with density = 1.0. One
of these sections should exist for every species we intend to use.

There are various types of profiles allowed, which are specified in the [species](Species.md) portion of the input deck. The currently available types are [*"standard"/"profile"*](#standard-profile), [*"constq"*](#constant-charge), [*"beamfocus"*](#beam-focus) and [*"file"*](#initialization-from-a-raw-particle-file), each of which are described below.

## Standard profile

This profile type is selected by setting `init_type` to *"standard"* or *"profile"* in the [species](Species.md) portion of the input deck. It accepts the following data:

- **density**, real, default, = 1.0
- **den_min**, real, default = 0.0
- **profile_type**(x_dim), character(\*), default = "uniform" (if num_x
  \> 0 then defaults to "piecewise-linear")

<!-- -->

- **num_x**, integer, default = -1
- **x**(num_x, x_dim), real, default = -Inf
- **fx**(num_x, x_dim), real, default = 0.0

<!-- -->

- **gauss_center**(x_dim), real, default = 0.0
- **gauss_sigma**(x_dim), real, default = Inf
- **gauss_range**(2,x_dim), real, default = {-Inf, Inf}
- **gauss_n_sigma**(x_dim), integer, default = 0

<!-- -->

- **channel_dir**, integer, default = 1
- **channel_r0**, float, default = 0.0
- **channel_depth**, float, default = 0.0
- **channel_size**, float, default = 0.0
- **channel_center**(x_dim-1), float, default = 0.0
- **channel_wall**, float, default = 0.0
- **channel_pos**(2), float, default = 0.0
- **channel_bottom**, float, default = 1.0

<!-- -->

- **sphere_center**(x_dim), float, default = 0.0
- **sphere_radius**, float, default = 0.0

<!-- -->

- **math_func_expr**, character(\*), default = "NO_FUNCTION_SUPPLIED!"

**density** - specifies a global multiplication factor for the density
profile. Regardless of which profile type you choose, the final density
value will be `density * profile` value (see below). This parameter can be used to set
the density for the sphere and uniform density profile types as these
default to 1.0 in the appropriate regions.

**den_min** - specifies the minimum density for injecting particles.
Particles are only injected when the specified density is above this threshold. This parameter is ignored if used in a
species connected to a `neutral` or `neutral_mov_ions` particle source.

**profile_type** - specifies the profile type to use.
There are two groups of density profiles that can be selected: (1) functions
that are separable into functions that depend only on one of the
coordinates, i.e., $f(x_1,x_2) = f_1(x_1) \times f_2(x_2)$ and (2)
arbitrary functions. If you choose the use the
first group, you can use piecewise-linear and Gaussian functions for any
of the $f_i(x_i)$. If you choose the second group, you can use one of the
following function types: uniform, channel, sphere or a specified
mathematical function. You cannot mix functions from the two groups.

If using separable functions, each component of the `profile_type` parameter must be set to
one of the following values:

- *"piecewise-linear"* - uses a piecewise-linear function defined by the
  `num_x`, `x` and `fx` parameters. See these parameters for details.
- *"gaussian"* - uses a Gaussian function defined by the `x_cent`, `f_cent`,
  `x_sigm`, and `x_range` parameters. See these parameters for details.

Here's an example for a 3D run using piecewise-linear functions in the
x1 and x2 direction and a Gaussian function in the x3 direction:

```text
  profile_type(1:3) = "piecewise-linear", "piecewise-linear", "gaussian",
```

If using an arbitrary function, you should set the first or all the
components of profile type to one of the following values:

- *"uniform"* - use a uniform density profile. The density value is set
  using the `density` parameter.
- *"channel"* - use a parabolic channel density profile. Channel
  parameters are set using the `channel_dir`, `channel_r0`,
  `channel_depth`, `channel_size`, `channel_center`, `channel_wall`,
  `channel_pos` and `channel_bottom`. See these parameters for details.
- *"sphere"* - use a spheric density profile. The density value is set
  using the `density` parameter. Sphere parameters are set using the
  `sphere_center` and `sphere_radius` parameters. See these parameters for
  details.
- *"math func"* - use the mathematical function supplied in the
  `math_func_expr` parameter to define the density profile. See this
  parameter for details.

Here's an example using a spheric density profile:

```text
  profile_type = "sphere",
```

**num_x**, **x**, **fx** - specify the parameters for the piecewise-linear
function we intend to use. The `num_x` parameter represents the number of points to
use in the piecewise-linear function (note that by setting this
parameter to a value \> 0 you also default your profile type to
piecewise-linear). The `x(k,i)` and `fx(k,i)` parameters represent the position and
density of the i-th direction for the k-th position. If `profile_type(i)`
is not equal to "piecewise-linear", the values specified are silently
ignored. The index `k` ranges from 1 to `num_x`,
and this range is the same for every `i` (ranging from 1 to `x_dim`). Here's
an example for a 6-point piecewise-linear function in the x2 direction:

```text
  x(1:6,2) = 0., 2.0, 2.1, 4.9, 5.0, 1000.,
  fx(1:6,2) = 0., 0., 1., 1., 0., 0.,
```

**gauss_center**, **gauss_sigma**, **gauss_range**, **gauss_n_sigma** - specify the
parameters for the Gaussian function we intend to use. The `gauss_center(i)`,
`gauss_sigma(i)` and `gauss_range(\*,i)` parameters represent the center, sigma and
range for the Gaussian on the i-th direction, respectively. If `profile_type(i)` is not
equal to *"gaussian"* the values specified are silently ignored. The
`gauss_range(k,i)` represent the lower (`k=1`) and upper (`k=2`) limits for
which the Gaussian function is evaluated. Outside these limits a value
of 0 is returned. The `gauss_n_sigma(i)` parameter, if specified, takes precedence over the `gauss_range(i)` parameter and sets the range to include `gauss_n_sigma(i)` number of sigmas on either side of the center. Here's an example for a Gaussian function in the x1
direction with sigma = 1.4142, centered about x1 = 10.0, with a peak
value of 2.0, defined for $5.0 < x_1 < 15.0$:

```text
  gauss_center(1) = 10.0,
  gauss_sigma(1) = 1.4142,
  gauss_range(1:2,1) = 5.0, 15.0,
  density = 2.0,
```

**channel_dir**, **channel_r0**, **channel_depth**,
**channel_size**, **channel_center**, **channel_wall**, **channel_pos**, **channel_bottom** -
specify the parameters for a parabolic channel. The `channel_dir` parameter specifies
channel symmetry axis. If `profile_type(1)` is not set to *"channel"*, the
values specified are silently ignored. The density of the channel will
be `channel_bottom + channel_depth * (||x_perp -
channel_center||/channel_r0)^2` for `||x_perp - channel_center|| <
channel_size/2`, where `x_perp` is the position perpendicular to the
symmetry axis. For `||x_perp - channel_center|| >= channel_size/2`,
the density will either be constant (`channel_wall <= 0.0`) or fall off
linearly to 0.0 along a distance of `channel_wall`. The `channel_pos(k)` parameter specifies
the beggining (`k=1`) and end (`k=2`) of the channel along the symmetry
axis. The `channel_center` parameter represents the coordinates perpendicular to the
symmetry axis, i.e., in 3D, if `channel_dir = 2`, then `channel_center = {
center x1, center x3 }`. Here's an example for a 2D run.

```text
  channel_dir = 1,
  channel_r0 = 0.349308,
  channel_depth = 32.727,
  channel_size = 2.09585,
  channel_center = 1.74654,
  channel_wall = 0.174654,
  channel_pos(1:2) = 15.0659, 10000.0,
  channel_bottom = 1.0,
```

**sphere_center**, **sphere_radius** - specify the parameters for a
spheric density profile. If `profile_type(1)` is not set to "sphere", the
values specified are silently ignored. The density value will be 1.0 for
`|| x - sphere_center || <= sphere_radius` and zero otherwise. The
`density` parameter should be used to change the density value inside the
sphere.

**math_func_expr** - specifies the mathematical function to be used to
define the profile. If `profile_type(1)` is not set to "math func", this
value is silently ignored. This expression can be a function of any of
the following, which represent the physical coordinates of the position
being calculated: `x1` (in 1D, 2D and 3D), `x2` (in 2D and 3D) and `x3` (in
3D). See the documentation on the analytical function parser for details
on the mathematical expression. Here's an example defining a pac-man-shaped density in 2D:

```text
  math_func_expr = "if((x1-6.4)^2+(x2-6.4)^2<6.4^2,
                       if( (abs(x2-6.4) < 4.5-x1),0.,1.),
                       0.0)",
```

## Constant charge

This profile type is selected by setting `init_type` to *"constq"* in the [species](Species.md) portion of the input deck. This profile should be used when particles are desired to have equal charge (rather than a fixed number of particles per cell). This profile only works for separable density functions and is thus incompatible with the *"channel"* and *"sphere"* profile types. It accepts all of the same parameters as the standard profile, except for the following additional parameter:

- **sample_rate**(x_dim), integer, default = 32

**sample_rate** - specify the number of sampling points per cell in each dimension.

## Beam focus

This profile type is selected by setting `init_type` to *"beamfocus"* in the [species](Species.md) portion of the input deck. This profile should be used when injecting a beam-like particle distribution with a specified Gaussian width, focus, etc. It accepts the following data:

- **density**, real, default, = 1.0
- **den_min**, real, default = 0.0

<!-- -->

- **gauss_center**(x_dim), real, default = 0.0
- **gauss_sigma**(x_dim), real, default = Inf
- **gauss_range**(2,x_dim), real, default = {-Inf, Inf}
- **gauss_n_sigma**(x_dim), integer, default = 0

<!-- -->

- **focal_dist**(3), real, default = 0.0
- **alpha**(3), real, default = 0.0
- **uth**(3), real, default = 0.0
- **gamma**, real, default = 1.0

**density** - specifies a global multiplication factor for the density
profile. Regardless of which profile type you choose the final density
value will be `density * profile` value (see below). This parameter can used to set
the density for the sphere and uniform density profile types as these
default to 1.0 in the appropriate regions.

**den_min** - specifies the minimum density for injecting particles.
Particles are only injected when the specified density (see profile
section) is above this threshold. This parameter is ignored if used in a
species connected to a neutral or neutral_mov_ions particle source.

**gauss_center**, **gauss_sigma**, **gauss_range**, **gauss_n_sigma** - specify the
parameters for the Gaussian function we intend to use. The `gauss_center(i)`,
`gauss_sigma(i)` and `gauss_range(\*,i)` parameters represent the center, sigma and
range for the Gaussian on the i-th direction, respectively. If `profile_type(i)` is not
equal to *"gaussian"* the values specified are silently ignored. The
`gauss_range(k,i)` represent the lower (`k=1`) and upper (`k=2`) limits for
which the Gaussian function is evaluated. Outside these limits a value
of 0 is returned. The `gauss_n_sigma(i)` parameter, if specified, takes precedence over the `gauss_range(i)` parameter and sets the range to include `gauss_n_sigma(i)` number of sigmas on either side of the center.

**focal_dist** - specify the focal plane distance after acceleration for each dimension.

**alpha** - specify the amount to rotate in the momentum-position phasespace along each dimension/component. The momentum is transformed according to `p(i) = p(i) - x(i) * alpha(i)` for the i-th dimension.

**uth** - specify the beam thermal momentum at the focal plane in each dimension.

**gamma** - specifies the Lorentz factor, $\gamma$, associated with the longitudinal momentum (`p1`) of the particles after acceleration.

## Initialization from a RAW particle file

This profile type is selected by setting `init_type` to *"file"* in the [species](Species.md) portion of the input deck. This profile initializes a particle distribution as specified by a RAW particle file. It accepts the following data:

- **file_name**, character(\*), default = ""
- **x1_offset**, real, default = 0.0
- **q_mult**, real, default = 1.0

**file_name** - specify the path to the file that contains the RAW particle data to be read in.

**x1_offset** - specify a position offset in the x1-direction for all particles. When specified, the position of a given particle will be set as `x1_part = x1_raw - x1_offset`.

**q_mult** - specify a charge scaling factor that multiplies the RAW particle charge at initialization.
