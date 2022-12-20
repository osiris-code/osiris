---
layout: single
classes: wide
title: Smooth
permalink: /reference/smooth
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the smoothing for electrical currents and
electromagnetic fields. If not present the code will not do any
smoothing. It accepts the following data:

- **type**(x_dim), character(\*), default = "none"
- **order**(x_dim), integer, default = 1
- **swfj**(3, max(order), x_dim), integer, default = 1
- **laserkdx_omega0**, real, default = -1.0
- **laserkdx_dir**, integer, default = 1

<!-- -->

- **digital_fc**, real, default = -∞
- **digital_A**, real, default = -∞

**type** specifies the type of smoothing to perform in each direction.
Valid values are:

- *"none*" - No smoothing. This is the same as not specifying this
  section.
- *"custom"* - Do a user-defined smooth using the **order**, and
  **swfj** parameters
- *"binomial"* - Do a binomial smoothing with a level set by the
  **order** parameter, i.e., apply a 1,2,1 stencil **order** times.
- *"compensated"* - Do a compensated binomial smoothing i.e. apply a
  1,2,1 stencil **(order-1)** times, followed by an appropriate
  compensator that eliminates $k^2$ dependency of the transfer
  function near *k* = 0. Using this filtering with **order = 2** is the
  same as a 4th-order tri-diagonal filter \[1,3\].
- *"5pass"* - Do a compensated 5-pass binomial smoothing i.e. apply a
  1,2,1 stencil 4 times, followed by a -5,14,-5 stencil.
- *"tri6"* - Do a 6th order tridiagonal filter \[3\].
- *"tri8"* - Do an 8th order tridiagonal filter \[3\].
- *"laserkdx"* - Do a binomial smoothing, i.e. apply a 1,2,1 stencil,
  followed by a compensator that does not attenuate the laser frequency
  for the given grid \[2\]. The compensator is similar to -1, 6, -1,
  where the value 6 is corrected according to the **laserkdx_omega0**
  parameter and the cell size in the smooth direction. Filtering is applied
  ONLY in the direction specified by 'laserkdx_dir' and is disabled in
  the other directions.
- *"digital"* - Tunable low pass filter \[4\]. This filter uses the
  'best' approximation to an ideal low pass filter. Filter behavior is
  controlled by the 'order', 'digital_fc' and 'digital_A' parameters.
  See below for details.

**order** specifies the number of passes for the 'binomial',
'compensated' and 'custom' filters and the number of terms for the
'digital' filter. The final filter stencil will use $2 \times order + 1$
points.

**swfj** specify the smoothing parameters for "custom" smooth. Smoothing
is done as a sequence of convolutions of 3-point stencils, that can be
defined independently for all directions. The number of stencils to be
applied for direction 'x_dim' is controlled by the *order('x_dim')*
parameter. `swfj(i,j,k)` with `i = {1,..,3}`, `j = {1, ..., order(x_dim)}`, and `k
= {1, ..., x_dim}` defines the stencil value for point i, to be used
in the smoothing pass j, for direction k. It is not necessary to
normalize the stencil coefficients; the code will do so automatically.

**laserkdx_omega0, laserkdx_dir** control the *laserkdx filter*.
**laserkdx_omega0** sets the laser frequency used in this smooth type.
This is the frequency that will not be attenuated by the smoothing.
**laserkdx_dir** sets the direction in which to do the filtering.

**digital_fc, digital_A** combined with **order** control the *digital*
filter. **digital_fc** specifies the cutoff frequency in units of the
Nyquist frequency, with values in the range \]0,1\[. The filter response
at this frequency will be 1/2. **digital_A** specifies the size of the
[Gibbs phenomenon](http://cnx.org/content/m10092/latest/) wiggles in
-db. 100 is a good choice. Increasing the value of **digital_A** lowers
the oscillations of the filter response close to the dropoff point but
also decreases the sharpness of the filter. The **order** parameter
controls the order of the filter. The resulting kernel will use
$2 \times order + 1$ points. Increasing the **order** parameter improves
the sharpness of the filter.

Here's a simple 2D example using binomial compensated smoothing in all
directions:

```text
smooth 
{
  order(1:2) = 2,
  type(1:2) = "compensated",
}
```

Here's a 2D/3D example where j is only smoothed in the x2 direction. j
is first smoothed in the x2 direction using window weights of (1,2,1)
then it is smoothed with a (2,0,2) window and finally it is smoothed
again with window weights of (1,1,1):

```text
smooth
{
  order(2) = 3,   
  
  !     -- all 3 points of window
  !     |   -- smoothing pass
  !     |   |  -- x2 direction
  !     |   |  |
  !     v   v  v
  swfj(1:3, 1, 2) = 1,2,1,   ! window weights for pass #1
  swfj(1:3, 2, 2) = 2,0,2,   ! window weights for pass #2
  swfj(1:3, 3, 2) = 1,1,1,   ! window weights for pass #3  
}
```

## Smoothing and Parallel Decomposition

The code will do the smoothing in a single
pass, by first calculating the larger smoothing window that corresponds
to the multi-pass smoothing specified, and then applying this larger
window only once. This is done mainly to avoid having to update the
current/fields on parallel jobs at the end of each smoothing pass which
results in better performance and allows the current smooth and boundary
update routines to be completely independent. However, this has the
downside of requiring extra guard cells on each node to accommodate the
larger smoothing window. Since Osiris only allows for guard cells
corresponding to a single remote node, this imposes a limit on the
minimal number of cells per node you can have in a given direction. With
the example above, the minimal number of cells per node along x2 is now
8, instead of the value whithout smoothing which is 1. There is no
workaround (yet) for this issue, and the code will not launch if you
specify an insufficient amount of cells in a given direction.

## References

\[1\] Birdsall & Langdon, "Plasma Physics via Computer Simulation",
appendix C, pp. 437-441, "Digital Filtering in One and Two Dimensions"

\[2\] S. Martins et al., "Numerical simulations of laser wakefield
accelerators in optimal Lorentz frames,” Comput Phys Commun, vol. 181,
no. 5, pp. 869–875, 2010

\[3\] J. Shang, "High-order compact-difference schemes for
time-dependent Maxwell equations", Journal of Computational Physics,
vol. 153, no. 2, pp. 312–333, 1999.

\[4\] R. Walraven, "Digital Filters", Proceedings of the Digital
Equipment User's Society, Fall, 1984.
