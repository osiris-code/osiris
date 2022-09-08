# Species

This section configures particle species settings and must be present in
the input file. One of these sections must exist for every species we
intend to use, either as a particle source for particles injected at the
begining of the simulations obeying some spatial profile and temperature
distribution, or as buffer to hold particles injected during the run
from some external source (e.g. cathode, neutral and neutral_mov_ions).
It accepts the following data:

- **name**, character(\*), default = "species \#id"
- **push_type**, character(\*), default = "simd" when available,
  "standard" otherwise
- **num_par_x**(x_dim), integer, default = 0
- **rqm**, real, default = 0.0

<!-- -->

- **num_par_max**, integer, default = 0

<!-- -->

- **n_sort**, integer, default = 25
- **push_start_time**, real, default = -1.0

<!-- -->

- **add_tag**, logical, default = .false.

<!-- -->

- **den_min**, real, default = 0.0
- **free_stream**, logical, default = .false.
- **subcycle**, bool, default = .false.

<!-- -->

- **q_real**, real, default = 0.0
- **if_collide**, logical, default = .false.

<!-- -->

- **if_like_collide**, logical, default = .false.

**name** specifies the name of this species to be used in diagnostic
output, both for file/directory names and data labels. If not specified
it will default to "species \#id", where \#id is the number of this
species (if this was the 3rd species to be specified, then \#id = 3). If
this species is connected to a cathode particle source then this
parameter defaults to "cathode \#id", where \#id is the number of the
cathode; if this species is connected to a neutral particle source then
this parameter defaults to "electrons \#id", where \#id is the number of
the neutral; and if this species is connected to a neutral_mov_ions
particle source then this parameter defaults to either "electrons \#id"
or "ions \#id" (depending on wether this is the electron or ion
species), where \#id is the number of the neutral_mov_ions particle
source. No two particles can have the same name (names that only differ
in case are considered equal. Space characters are considered equal to
underscore "_" characters).

**push_type** specifies the algorithm used for pushing the particles.
Currently available options are:

- "simd" - Push the particles using a SIMD accelerated pusher. This is
  the default when available. Some features (e.g. free streaming) are
  not available when using this pusher.
- "standard" - Push the particles using the standard Fortran pusher.
- "pgc" - Push the particles using the ponderomotive guiding center
  algorithm pusher.

**num_par_x** specifies the number of particles per cell to use in each
direction. The total number of particles per cell will be the product of
all the components of num_par_x

**rqm** specifies the 'reciprocal charge over mass ratio' (i.e. m/q) in
regard to electron mass and absolute charge for this species. Electrons
would be represented by rqm = -1.0 and protons would be represented by
choosing rqm = 1836.1527.

**num_par_max** specifies the particle buffer size for each node. When
not set the code will try to determine this value automatically. If for
some reason the number of particles on each node exceeds the supplied
value the code will attemp to resize the buffer to accomodate the extra
particles.

**n_sort** specifies the number of iteration between each particle sort.
This operation results in better cache use and significant performance
enhancement for long runs. The default value of 25 represents a good
balance between the extra time it takes to sort the particles and the
performance gain. Setting this parameter to 0 turns off particle sorting
for this species.

**push_start_time** specifies the time in simulation units at which to
start pushing the particles. The species will not be pushed unless the
simulation time is greater than this parameter.

**add_tag** specifies whether the code should add a unique tag to each
particle in this species so it can be identified for diagnostic
purposes.

**den_min** specifies the minimum density for injecting particles.
Particles are only injected when the specified density (see profile
section) is above this threshold. This parameter is ignored if used in a
species connected to a neutral or neutral_mov_ions particle source.

**free_stream** setting this parameter to .true. will prevent the code
from changing the species velocity, so particles will stream freely.
Electrical current is still deposited.

**subcycle** specifies whether to use time averaged fields to push this
species. See also the n_subcycle parameter in the el_mag_fld section.
Note that this is still at an experimental stage, and the physics hasn't
been properly tested. The code shoud, however, run ok

**q_real** specifies the charge of this species normalized to the
elementary charge. Electrons would be represented by q_real = -1.0 and
protons would be represented by q_real = 1.0. This parameter just needs
to be specified if this species makes Coulomb collisions.

**if_collide** setting this parameter to .true. will allow this species
to make binary Coulomb collisions.

**if_like_collide** setting this parameter to .true. will allow this
species to make binary Coulomb collisions with itself.

Here's an example of a species section for a 2D run. It will be named
"electrons", it uses electrons or a particle with the same q/m, 8
particles per cell, 4 in the x1 direction and 2 in the x2 direction, and
the buffer size will be set automatically.

```text
species
{
 name = "electrons",
 rqm=-1.0,
 num_par_x(1:2) = 4, 2,
}
```
