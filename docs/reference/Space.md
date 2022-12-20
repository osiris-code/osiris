---
layout: single
classes: wide
title: Space
permalink: /reference/space
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the spatial information and moving window
settings of the simulation and must be present in the input file. It
accepts the following data:

- **xmin**(x_dim), float, default = 0.0
- **xmax**(x_dim), float, default = 0.0
- **if_move**(x_dim) bool, default = .false.
- **move_u**, float, default = -1.0

**xmin**, **xmax** - specify the lower and upper boundaries of the global
simulation space at the beggining of the simulation in normalized units.
Note that for cylindrical coordinates, the lower boundary of the radial
dimension will always be set to -*dr*/2, where *dr* is the radial cell
size, overriding user input.

**if_move** - specifies a whether the code should use a moving window at
the speed of light in the specified directions. Setting this to .true.
will override any boundary conditions for that direction you specify for
Electro-Magnetic fields or particle species, except for the periodic
boundaries set on the node_conf section.

**move_u** - specifies the normalized momentum of the moving window (i.e. beta * gamma).
Note: in principle, negative momenta could be allowed, but are not currently implemented.

Here's an example of a space section for a 2D run that sets the global
simulation boundaries to \[0.0,15.0\] in the x1 direction and \[0.0,
3.5\] in the x2 direction. The simulation will use a moving window along
the x1 direction.

```text
space
{
  xmin(1:2) =  0.0, 0.0,
  xmax(1:2) =   15.0, 3.5,
  if_move(1:2) = .true.,.false.,
}
```
