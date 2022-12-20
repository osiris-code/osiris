---
layout: single
classes: wide
title: Collisions
permalink: /reference/collisions
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures binary Coulomb collisions settings.

Please note that the OSIRIS collision package is still in the
'development' stage (where it possibly will be forever.) First, as
described below, most options besides the default ones are currently
undertested if not nonfunctional. Anyone wishing to use those functions
is advised to work with the development team to finish implementing and
benchmarking them.

Furthermore, in a broader sense, the correct and efficient addition of
collisions to the PIC algorithm remains an open problem. To a certain
extent doing so inherently 'breaks' the PIC algorithm, in that it
inserts a somewhat arbitrary force into the particle-field push loop.
Anyone using this package is advised to familiarize themselves with the
basic functioning of the code (an effort has been made to include
elucidating inline documentation); a partial bibliography has also been
included in the header of the source file. That said, if the user sets
the default options, and sufficiently resolves the collision physics
(both in terms of collision steps per collision period and colliding
particles per collision cell), the package should 'work' correctly.

The input deck can include the following data:

- **n_collide**, integer, default = 0
- **nx_collision_cells**(x_dim), integer, default = 1
- **coulomb_logarithm_automatic**, logical, default = .true.
- **coulomb_logarithm_value**, real, default = 0.0
- **collision_model**, string, default = 'Sentoku'
- **cross_section_correction**, string, default = 'None'
- **timestep_cross_corr**, logical, default = .false.
- **frame_correction**, string, default = 'None'
- **root_finder**, string, default = 'None'

Note it is also necessary to set either the **n0** or the **omega_p0**
parameters in the [simulation section](simulation) to specify the
reference density; and the **q_real**, **if_collide**, and
**if_like_collide** parameters in [species section](species) for any species which are
colliding.

**n_collide** specifies the number of iterations between collisions. If
set to 0 collisions will be turned off. This parameter must be a
multiple of **n_sort** (see the [species section](species)).

**nx_collision_cells** specifies the number of simulation cells per
collisional cell in each direction. The PIC grid must be a multiple of
the collisional grid on each node (i.e. the total grid divided by the
number of nodes in a given direction.) Currently the load balancer does
not respect the collision cells, so nx_collision_cells should be set to
1 in any balanced direction.

**coulomb_logarithm_automatic** specifies whether to self-consistently
calculate the Coulomb logarithm from the local plasma density and
temperature (default) or use a fixed Coulomb logarithm value.

**coulomb_logarithm_value** specifies the Coulomb logarithm value to be
used if coulomb_logarithm_automatic is set to false.

------------------------------------------------------------------------

The following options preserve all of the functionality that was present
as pre-compiler options in the old code (pre 2015). However the options
other than the defaults have not been tested thoroughly and caution is
urged in using any of them. Most users will want to ignore these fields.
All strings are case-insensitive.

**collision_model** specifies which collision model to use. Accepted
options are 'Sentoku', 'Nanbu', 'Nanbu_rel', and 'Isotropic'. Note that
'Nanbu' and 'Nanbu_rel' are only rudimentarily implemented at the moment
(see **root_finder** below) and so if it is desired to use them that
section of the code will need to be developed. Also note that
'Isotropic' scatters into 4π at every collision step and so is only
meant for testing.

**cross_section_correction** specifies the method to transform the
effective cross section in a way which maintains Lorentz invariance (see
Peano et al. 2009). Accepted options are 'None', 'Rejection', and
'Scattering.'

**timestep_cross_corr** is to be used when 'cross_section_correction' is
set to anything other than 'None;' it multiplies the effective cross
section by 2, to correct for the reduction of cross section in the lab
frame caused by the cross section correction algorithm.

**frame_correction** accepts the strings 'none', 'time', and
'time_density'. This option increases the cross section by the
relativistic gamma factor of the more massive particle as seen in the
lab frame, either by one factor ('time'), or by two factors
('time_density'). This is preserved from the earlier version of the code
but its use is not well documented. Unless you know you need this, it is
recommended to leave this as 'none'.

**root_finder** specifies which algorithm to use to invert the Nanbu
equation for the cross section. It accepts 'Newton' for a Newton solver;
'Linear' for a linear-interpolation look-up table; 'Cubic' for a cubic
interpolation look-up table; and 'None'. However, if using either Nanbu
collision model, all options other than 'Newton' will result in a
run-time error; 'Linear' and 'Cubic' are not currently implemented (the
look-up tables need to be written) and 'None' is meant as an
inconsistent default option to make sure the user actively chooses the
algorithm. In addition, the Newton solver is very inefficient and is not
advised for production runs. Therefore if the user wishes to use a Nanbu
collision model, it is highly recommended that a more efficient
algorithm be implemented (either writing the look-up tables for the
interpolating solvers, or implementing the improved method of Perez et
al. 2012)

------------------------------------------------------------------------

Here is an example of a collisions section.

```text
collisions
{
 n_collide = 1,  ! collide every time step
 nx_collision_cells(1:2) = 8,8, ! each collision cells has 8x8 PIC cells
   
 coulomb_logarithm_automatic = .false., ! fixed Coulomb logarithm
 coulomb_logarithm_value = 4.0,
}
```
