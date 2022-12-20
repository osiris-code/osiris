---
layout: single
classes: wide
title: Species
permalink: /reference/species
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures particle species settings and must be present in
the input file. One of these sections must exist for every species we
intend to use, either as a particle source for particles injected at the
begining of the simulations obeying some spatial profile and temperature
distribution, or as buffer to hold particles injected during the run
from some external source (e.g., [cathode](Cathode.md), [neutral](Neutrals.md) and [neutral_mov_ions](Neutrals_with_Moving_Ions.md)).
It accepts the following data:

- **name**, character(\*), default = "species \#id"
- **push_type**, character(\*), default = "simd" when available,
  "standard" otherwise
- **rqm**, real, default = 0.0

<!-- -->

- **num_par_x**(x_dim), integer, default = -1
- **num_par_max**, integer, default = 0
- **tot_par_x**(x_dim), integer, default = -1

<!-- -->

- **n_sort**, integer, default = 25
- **push_start_time**, real, default = -1.0
- **add_tag**, logical, default = .false.
- **free_stream**, logical, default = .false.

<!-- -->

- **q_real**, real, default = 0.0
- **if_collide**, logical, default = .false.
- **if_like_collide**, logical, default = .false.

<!-- -->

- **num_pistons**, integer, default = 0

<!-- -->

- **init_fields**, logical, default = .false.
- **init_type**, character(\*), default = "standard"

<!-- -->

- **iter_tol**, real, default = 1.0e-3
- **rad_react**, logical, default = .false.

**name** - specifies the name of this species to be used in diagnostic
output, both for file/directory names and data labels. If not specified
it will default to "species \#id", where \#id is the number of this
species (if this was the third species specified, then \#id = 3). If
this species is connected to a cathode particle source, then this
parameter defaults to "cathode \#id", where \#id is the number of the
cathode; if this species is connected to a neutral particle source then
this parameter defaults to "electrons \#id", where \#id is the number of
the neutral; and if this species is connected to a neutral_mov_ions
particle source then this parameter defaults to either "electrons \#id"
or "ions \#id" (depending on whether this is the electron or ion
species), where \#id is the number of the neutral_mov_ions particle
source. No two species can have the same name (names that only differ
in case are considered equal, and space characters are considered equal to
underscore "\_" characters).

**push_type** - specifies the algorithm used for pushing the particles.
Currently available options are:

- *"simd"* - Push the particles using a SIMD accelerated pusher that
  employs the standard Boris algorithm \[1\]. This is
  the default when available. Some features (e.g. free streaming) are
  not available when using this pusher.
- *"standard"* - Push the particles using the standard Fortran pusher
  using the Boris algorithm \[1\].
- *"vay"* - Push the particles using the Vay algorithm \[2\], which eliminates a spurious force on relativistic particles when using the Boris algorithm.
- *"cond_vay"* - A conditional Vay pusher, where particles with $\gamma \leq 5$ are pushed using the Boris algorithm, and otherwise with the Vay algorithm (the $\gamma$ condition is somewhat arbitrary and currently hard-coded into `dudt_cond_vay`).  The Boris algorithm is usually more accurate for low-energy particles.
- *"fullrot"* - This algorithm calculates the magnetic field rotation of the Boris algorithm exactly (rather than second-order accurate).  However, the electric field half-pushes remain only first-order accurate in the Boris algorithm \[1,3\].
- *"euler"* - This algorithm implements a Boris-type pusher where the magnetic rotation is implemented using the [Euler-Rodrigues formula](https://en.wikipedia.org/wiki/Euler–Rodrigues_formula), which gives exact results for arbitrarily high magnetic field values.
- *"cary"* - Push the particles using the algorithm proposed by Higuera and Cary \[4\].  This method is volume-preserving (like Boris) and gives the correct $E \times B$ force (like Vay).
- *"exact"* or *"analytic"* - Performs an exact integration of the particle momenta assuming constant fields for a given time step \[5\].  The `rad_react` parameter can be set to `.true.` to include a perturbation from radiation reaction with this pusher, and the `iter_tol` parameter is used to specify the root-finding convergence threshold.
- *"exact-rr"* or *"analytic-rr"* - Performs an exact integration of the particle momenta based on the semi-classical form of radiation reaction in the Landau-Lifshitz equation, assuming constant fields for a given time step \[5\]. The `n0` or `omega_p0` parameter should be specified in the [simulation](Simulation.md) portion of the input deck, and the `q_real` parameter should be set within this section. For the strong-field regime, the *“exact-rr”* pusher is recommended due its improved accuracy compared to the *“exact”* pusher with `rad_react` set to `.true.`.
- *"radcool"* - This algorithm takes into account classical radiation reaction using the Landau-Lifshitz equation of motion, where particle motion is integrated with the fourth-order Runge-Kutta method (we neglect the explicit temporal derivatives in the L-L equation) \[6\]. If specified, the `n0` or `omega_p0` parameter should be specified in the [simulation](Simulation.md) portion of the input deck.

**rqm** - specifies the 'reciprocal charge over mass ratio' (i.e., m/q) in
regard to electron mass and absolute charge for this species. Electrons
would be represented by `rqm = -1.0` and protons would be represented by
choosing rqm = 1836.1527.

**num_par_x** - specifies the number of particles per cell to use in each
direction. The total number of particles per cell will be the product of
all the components of `num_par_x`.

**num_par_max** - specifies the particle buffer size for each node. When
not set, the code will try to determine this value automatically. If for
some reason the number of particles on each node exceeds the supplied
value, the code will attempt to resize the buffer to accommodate the extra
particles.

**tot_par_max** - specifies the total number of particles in each dimension of the simulation (across all grids). This quantity will be overridden if `num_par_x` is specified.

**n_sort** - specifies the number of iterations between each particle sort.
This operation results in better cache use and significant performance
enhancement for long runs. The default value of 25 represents a good
balance between the extra time it takes to sort the particles and the
performance gain. Setting this parameter to 0 turns off particle sorting
for this species.

**push_start_time** - specifies the time in simulation units at which to
start pushing the particles. The species will not be pushed unless the
simulation time is greater than this parameter.

**add_tag** - specifies whether the code should add a unique tag to each
particle in this species so it can be identified for diagnostic
purposes (e.g., for particle tracking).

**free_stream** - setting this parameter to `.true.` will prevent the code
from changing the species velocity, so particles will stream freely.
Electrical current is still deposited.

**q_real** - specifies the charge of this species normalized to the
elementary charge. Electrons would be represented by `q_real = -1.0`, and
protons would be represented by `q_real = 1.0`. This parameter only needs
to be specified if this species makes Coulomb collisions or is using the exact pusher with radiation reaction.

**if_collide** - setting this parameter to `.true.` will allow this species
to make binary Coulomb collisions. If specified, the `q_real` parameter should be set within this section.

**if_like_collide** - setting this parameter to `.true.` will allow this
species to make binary Coulomb collisions with itself. If specified, the `q_real` parameter should be set within this section.

**num_pistons** - specifies the number of pistons to be read in.

**init_fields** - setting this parameter to `.true.` will initialize the fields associated with the initial density/momentum particle distribution using an electrostatic field solver.

**init_type** - specifies the type of profile initialization for this species.  Currently available options are:

- *"standard"/"profile"* - Standard initialization from a specified density profile.
- *"constq"* - Injects particles into an area defined by certain grid cell indices using a fixed charge per particle (as opposed to a fixed number of particles per cell).
- *"beamfocus"* - Used to inject a particle beam with specified Gaussian width, focus, etc.
- *"file"* - Used to inject particles from a specified RAW particle file.

**iter_tol** - specifies the convergence threshold of root-finding routines in the exact pusher.

**rad_react** - setting this parameter to `.true.` will include radiation reaction by correcting the particle momenta each time step. This parameter should only be set when *“exact”* particle pusher is used. If `rad_react` is specified, the `n0` or `omega_p0` parameter should be specified in the [simulation](Simulation.md) portion of the input deck, and the `q_real` parameter should be set within this section.

Here's an example of a species section for a 2D run. It will be named
"electrons", it uses electrons or a particle with the same q/m, 8
particles per cell, 4 in the x1 direction and 2 in the x2 direction, and
the buffer size will be set automatically.

```text
species
{
  name = "electrons",
  rqm = -1.0,
  num_par_x(1:2) = 4, 2,
}
```

## References

\[1\] J. P. Boris, "Relativistic plasma simulation—Optimization of a hybrid code," 4th Conference on Numerical Simulation of Plasmas, pp. 3–67, 1970.

\[2\] [J.-L. Vay, "Simulation of beams or plasmas crossing at relativistic velocity," Phys Plasmas, vol. 15, no. 5, p. 056701, 2008.](https://doi.org/10.1063/1.2837054)

\[3\] [V. K. Decyk, et al., "An Analytic Boris Pusher for Plasma Simulation," arXiv preprint arXiv:2204.13563, 2022.](https://doi.org/10.48550/arXiv.2204.13563)

\[4\] [A. V. Higuera and J. R. Cary, "Structure-preserving second-order integration of relativistic charged particle trajectories in electromagnetic fields," Phys Plasmas, vol. 24, no. 5, p. 052104, 2017.](https://doi.org/10.1063/1.4979989)

\[5\] [F. Li, et al., "Accurately simulating nine-dimensional phase space of relativistic particles in strong fields," J Comput Phys, vol. 438, p. 110367, 2021.](https://doi.org/10.1016/j.jcp.2021.110367)

\[6\] [Vranic, M., et al., “Classical radiation reaction in particle-in-cell simulations,” Comput Phys Commun, vol. 204, pp. 141–151, 2016.](https://doi.org/10.1016/j.cpc.2016.04.002)
