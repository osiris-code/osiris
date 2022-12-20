---
layout: single
classes: wide
title: Simulation
permalink: /reference/simulation
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures global simulation parameters and is optional. It
accepts the following data:

- **random_seed**, integer, default = 0
- **random_algorithm**, string, default = "default"
- **omega_p0**, real, default = 0.0
- **n0**, real, default = 0.0
- **gamma**, real, default = 1.0
- **ndump_prof**, integer, default = 0
- **wall_clock_limit**, string, default = ""
- **wall_clock_check**, integer, default = -1
- **wall_clock_checkpoint**, logical, default = .true.
- **prof_start**, integer, default = -1
- **prof_end**, integer, default = huge(prof_start)
- **file_format**, string, default = "hdf5"
- **parallel_io**, string, default = "mpi"
- **enforce_rng_constancy**, logical, default = .false.

**random_seed** specifies the seed used by the random number generator.
If not set it defaults to 0. If set the random number generator seed
used by OSIRIS in each parallel node will be random_seed + node_id (or
just random_seed in serial runs).

**random_algorithm** algorithm to use for global random number generation.
Currently available options are:

- *"mt"* - Fortran implementation of the [Mersenne Twister random number generator](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html).

- *"r250"* - F03 Implementation of the R250 Pseudo-random number generator \[1\], \[2\]

- *"cmwc"* - Complimentary Multiply With Carry Random Number Generator \[3\]

- *"mwc"* - Marsaglia MWC (multiply-with-carry) Pseudo Random Number Generator

- *"kiss"* - KISS Pseudo Random Number Generator \[3\]

- *"hash"* - Implements the Hashing RNG from \[4\]

- *"default"* - Uses Mersenne Twister, above.

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
in a format "hours:minutes:seconds". When this wall clock time is reached the simulation
will shut down. See also 'wall_clock_check'.

**wall_clock_check** specifies the frequency at which to check if wall
time limit has been reached. See also 'wall_clock_limit'.

**wall_clock_checkpoint** specifies if the code should dump
checkpointing information when stopping due to wall clock limit.

**prof_start** specifies the time step at which profiling dumps should begin.

**prof_end** specifies the time step at which profiling dumps should end.

**file_format** specifies the file format to be used for simulation output.
Currently available options are:

- *"hdf5"* - Output files using the [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/). Using this option requires compiling OSIRIS with HDF5 support. When supported, this is the default.
- *"zdf"* - Output files using the [ZDF format](https://github.com/zpic-plasma/zpic/tree/main/zdf). This option generally leads to smaller output files and better performance but it is not wildly supported by visualization and data analysis tools.

**parallel_io** Specifies the type of parallel I/O used for simulation output.
Currently available options are:

- *"mpi"* (default) - Data is saved by root node and transferred from other nodes using MPI. Communication and file I/O are overlapped for better performance. Only root node needs access to the file system.
- *"indep"* (not supported by HDF5 files) - Data is saved by all nodes independently. Requires direct access to file system by all nodes.
- *"mpiio-indep"* Data is saved using MPI parallel I/O (MPI-IO) using independent writes. If using HDF5 this requires that the HDF5 library was compiled with parallel I/O support.
- *"mpiio-coll"* Data is saved using MPI parallel I/O (MPI-IO) using collective writes. If using HDF5 this requires that the HDF5 library was compiled with parallel I/O support.

**enforce_rng_constancy** If false, an identical iput deck will give different result
due to random number generation (RNG) depending on a given run's spatial decomposition
(e.g. different number/load of MPI NODES).
If true, an identical input deck identical results (within roundoff) no matter the run's MPI setup
(i.e. running a simulation on 1 node will give same results no matter how many MPI nodes used).
This is comes at a small perforance cost, but is usefull for testing, debugging the and subtraction-technique.

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

## References

\[1\] Kirkpatrick, S., and E. Stoll, A Very Fast Shift-Register Sequence Random Number Generator, Journal of Computational Physics (1981), v. 40. p. 517

\[2\] Maier, W.L., 1991; A Fast Pseudo Random Number Generator, Dr. Dobb's Journal, May, pp. 152 - 157

\[3\] Marsaglia, George (2003) "Random Number Generators," Journal of Modern Applied Statistical Methods: Vol. 2: Iss. 1, Article 2.

\[4\] Saul Teukolsky, William H. Press, and William T. Vetterling, "Numerical Recipes", 2007
