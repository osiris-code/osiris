# Profile

This section configures the density profile for the particle species and
must be present in the input file and is optional. If not specified the
code will default to an uniform density profile with density = 1.0. One
of these sections should exist for every species we intend to use. It
accepts the following data:

- **density**, real, default, = 1.0
- **profile_type**(x_dim), character(\*), default = "uniform" (if num_x
  \> 0 then defaults to "piecewise-linear")
- **num_x**, integer, default = -1
- **x**(num_x, x_dim), real, default = -Inf
- **fx**(num_x, x_dim), real, default = 0.0
- **gauss_center**(x_dim), real, default = 0.0
- **gauss_sigma**(x_dim), real, default = Inf
- **gauss_range**(2,x_dim), real, default = {-Inf, Inf}
- **channel_dir**, integer, default = 1
- **channel_r0**, float, default = 0.0
- **channel_depth**, float, default = 0.0
- **channel_size**, float, default = 0.0
- **channel_center**(x_dim-1), float, default = 0.0
- **channel_wall**, float, default = 0.0
- **channel_pos(2)**, float, default = 0.0
- **channel_center**, float, default = 0.0
- **channel_bottom**, float, default = 1.0
- **sphere_center**(x_dim), float, default = 0.0
- **sphere_radius**, float, default = 0.0
- **math_func_expr**, character(\*), default = "NO_FUNCTION_SUPPLIED!"

There are two groups of density profiles that can be selected: functions
that are separable into functions that depend only on one of the
coordinates i.e. $f(x_1,x_2) = f_1(x_1) \times f_2(x_2)$ belong the first group, and
arbitrary functions fall into the second. If you choose the use the
first group you can use piecewise linear and gaussian functions for any
of the $f_i(x_i)$. If you choose the second group you can use one of the
following function types: uniform, channel, sphere or a specified
mathematical function. You cannot mix functions from the two groups.

**density** specifies a global multiplication factor for the density
profile. Regardless of which profile type you choose the final density
value will be `density * profile` value. This parameter can used to set
the density for the sphere and uniform density profile types as these
default to 1.0 in the appropriate regions.

**profile_type** specifies the profile type to use. If using separable
functions each component of the profile_type parameter must be set to
one of the following values:

- *"piecewise-linear"* - uses a piecewise linear function defined by the
  num_x, x and fx parameters. See these parameters for details.
- *"gaussian"* - uses a gaussian function defined by the x_cent, f_cent,
  x_sigm, and x_range parameters. See these parameters for details.

Here's an example for a 3D run using piecewise linear functions in the
x1 and x2 direction and a gaussian function in the x3 direction:

```text
profile_type(1:3) = "piecewise-linear", "piecewise-linear", "gaussian", 
```

If using arbitrary function you should set the first or all the
components of profile type to one of the following values:

- *"uniform"* - use a uniform density profile. The density value is set
  using the density parameter.
- *"channel"* - use a parabolic channel density profile. Channel
  parameters are set using the channel_dir, channel_bottom, channel_r0,
  channel_depth, channel_size, channel_center, channel_wall, and
  channel_pos. See these parameters for details.
- *"sphere"* - use a spheric density profile. The density value is set
  using the density parameter. Sphere parameters are set using the
  sphere_center, and sphere_radius parameters. See these parameters for
  details.
- *"math func"* - use the mathematical function supplied in the
  math_func_expr parameter to define the density profile. See this
  parameter for details.

Here's an example using a spheric density profile:

```text
profile_type = "sphere"
```

**num_x**, **x**, **fx** specify the parameters for the piecewise linear
function we intend to use. **num_x** represents the number of points to
use in the piecewise linear function (Note that by setting this
parameter to a value \> 0 you also default your profile type to
piecewise linear). **x**(k,i) and **fx**(k,i) represent the position and
density on the i-th direction for the k-th position. If profile_type(i)
is not equal to "piecewise-linear" the values specified are silently
ignored. k ranges from 1 to num_x (which is set in the num_x section)
and this range is the same for every i (ranging from 1 to x_dim). Here's
an example for a 6 point piecewise linear function in the x2 direction:

```text
x(1:6,2) = 0., 2.0, 2.1, 4.9, 5.0, 1000.,
fx(1:6,2) = 0., 0., 1., 1., 0., 0.,
```

**gauss_center**, **gauss_sigma**, **gauss_range** specify the
parameters for the gaussian function we intend to use. gauss_center(i),
gauss_sigma(i) and gauss_range(\*,i) represent the center, sigma and
range for the gaussian on the i-th direction. If profile_type(i) is not
equal to "gaussian" the values specified are silently ignored. The
gauss_range(k,i) represent the lower (k=1) and upper (k=2) limits for
which the gaussian function is evaluated. Outside these limits a value
of 0 is returned. Here's an example for a gaussian function in the x1
direction with sigma = 1.4142, centered about x1 = 10.0, with a peak
value of 2.0, defined for $5.0 < x_1 < 15.0$:

```text
gauss_center(1) = 10.0,
gauss_sigma(1) = 1.4142,
gauss_range(1:2,1) = 5.0, 15.0,
density = 2.0,
```

**channel_dir**, **channel_bottom**, **channel_r0**, **channel_depth**,
**channel_size**, **channel_center**, **channel_wall**, **channel_pos**
specify the parameters for a parabolic channel. channel_dir specifies
channel symmetry axis. If profile_type(1) is not set to "channel" the
values specified are silently ignored. The density of the channel will
be channel_bottom + channel_depth \* (\|\|x_perp -
channel_center\|\|/channel_r0)^2 for \|\|x_perp - channel_center\|\| \<
channel_size/2, where x_perp is the position perpendicular to the
simmetry axis. For \|\|x_perp - channel_center\|\| \>= channel_size/2,
the density will either be constant (channel_wall \<= 0.0) or fall off
linearly to 0.0 along a distance of channel_wall. channel_pos(k) specify
the beggining (k=1) and end (k=2) of the channel along the symmetry
axis. channel_center represents the coordinates perpendicular to the
symmetry axis i.e. in 3D, if channel_dir = 2 then channel_center = {
center x1, center x3 }. Here's an example for a 2D run.

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

**sphere_center**, **sphere_radius** specify the parameters for a
spheric density profile. If profile_type(1) is not set to "sphere" the
values specified are silently ignored.The density value will be 1.0 for
\|\| x - sphere_center \|\| \<= sphere_radius and zero otherwise. The
density parameter should be used to change the density value inside the
sphere.

**math_func_expr** specifies the mathematical function to be used to
define the profile. If profile_type(1) is not set to "math func" this
value is silently ignored. This expression can be a function of any of
the following, which represent the physical coordinates of the position
being calculated: x1 (in 1D, 2D and 3D), x2 (in 2D and 3D) and x3 (in
3D). See the documentation on the analytical function parser for details
on the mathematical expression. Here's an example defining a pac-man
shaped density in 2D:

```text
math_func_expr = "if((x1-6.4)^2+(x2-6.4)^2<6.4^2,
                                   if( (abs(x2-6.4) < 4.5-x1),0.,1.),
                                  0.0)",
```
