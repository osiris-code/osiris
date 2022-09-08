# Simulation

This section configures global simulation parameters and is optional. It
accepts the following data:

- **algorithm**, string, default = "standard"
- **random_seed**, integer, default = 0
- **omega_p0**, real, default = 0.0
- **n0**, real, default = 0.0
- **gamma**, real, default = 1.0
- **ndump_prof**, integer, default = 0
- **wall_clock_limit**, string, default = ""
- **wall_clock_check**, integer, default = -1
- **wall_clock_checkpoint**, logical, default = .true.

**algorithm** specifies the algorithm to be used in the simulation. If
set to "pgc" the ponderomotive guiding center solver (PGC) algorithm and
corresponding particle pusher are used. PGC is currently implemented in
2D geometries only. Please read node configuration and space
configuration sections for additional details on PGC usage. If set to
"hd-hybrid" the hybrid algorithm is used.

**random_seed** specifies the seed used by the random number generator.
If not set it defaults to 0. If set the random number generator seed
used by OSIRIS in each parallel node will be random_seed + node_id (or
just random_seed in serial runs).

**omega_p0** specifies the reference simulation plasma frequency for the
simulation in units of \[rad/s\]. The user must set this value when
doing simulations that require real (not simulation) units, such as
simulations involving ionization or collisions. Setting a reference
density, **n0**, overrides this value. See also the **n0** option.

**n0** specifies the reference simulation plasma density for the
simulation in units of \[cm^-3\]. Setting this value overrides the
**omega_p0** option, calculating its value from the supplied value. See
also the **omega_p0** option.

**gamma** specifies the relativistic factor of the simulation frame.
This factor will transform all zpulses to the boosted frame, and will be
used to automatically calculate the 'laserkdx' smooth coefficients.

**ndump_prof** specifies the frequency for profiling (timing)
information dumps.

**wall_clock_limit** specifies a time limit to shut down the simulation,
in a format "h:m:s". When this wall clock time is reached the simulation
will shut down. See also 'wall_clock_check'.

**wall_clock_check** specifies the frequency at which to check if wall
time limit has been reached. See also 'wall_clock_limit'.

**wall_clock_checkpoint** specifies if the code should dump
checkpointing information when stopping due to wall clock limit.

## Example

Here's an example of a simulation section defining a random number seed
of 1237 and a reference density of 2.7e17 cm^-3.

```text
simulation 
{
  random_seed = 1237, 
  n0 = 2.7e17,  ! [cm^-3]
}
```
