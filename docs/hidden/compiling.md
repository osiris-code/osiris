# Compiling OSIRIS

## Introduction

Compilation of OSIRIS on any system is done in 2 steps:

1.  **Configuration** - Selecting the appropriate Fortran/C compilers
    and compiler flags, MPI and HDF5 libraries and compilation options,
    and compilation target.
2.  **Compilation** - Creates the OSIRIS binary for execution.

Almost all OSIRIS options can be selected in the input file (referred to
as runtime options), so usually you will only have to compile the code
once. Options that are only available when compiling the code (compile
time options) are for example the dimensionality of the code (1D/2D/3D)
and numerical precision (single/double).

## Pre-requisites

To compile OSIRIS the following must be available in the target system:

Compilers/Tools:

- A Fortran 9X compiler (for x86 systems Intel ifort is recommended).
  OpenMP support is optional (but recommended).
- A C compiler (for x86 systems Intel ifort is recommended)
- gcc - This is used as a preprocessor for the fortran code
- A gnu compatible Make (AIX make is not compatible, but gmake is
  usually available)

Libraries:

- MPI - Any implementation should do:
  [MPICH2](http://www.mcs.anl.gov/research/projects/mpich2/), [Open
  MPI](http://www.open-mpi.org) and
  [MVAPICH](http://mvapich.cse.ohio-state.edu) are known to work. If the
  implementations doesn't support MPI_IN_PLACE operations these can be
  disabled in source/os-config.h. MPE can also be used for logging and
  profiling if available.
- [HDF5](http://www.hdfgroup.org/HDF5/) . Version 1.8 is recommended,
  although 1.6 should also work. Support for parallel I/O is optional.
- [PAPI](http://icl.cs.utk.edu/papi/) (optional). PAPI can optionally be
  used for detailed code profiling. Only the C interface is required.
- [SIONlib](http://www2.fz-juelich.de/jsc/sionlib) (optional) - SIONlib
  can be used to reduce the number of files used for writing restart
  information.

## Configuring compilation

OSIRIS compilation is configured by issuing the "configure" command at
the root of the OSIRIS source code:

` $ ./configure [-h] [-s system] [-d dimension] [-t type]`

Available options are:

- -h - Display a help message and a list of available systems.
- -s - Set the system configuration. System configuration is defined in
  a file named osiris_sys.*system_name* saved in the ./config directory.
  See below for details
- -d - The dimensionality of the algorithm. Can be 1, 2 or 3 (defaults
  to 2).
- -t - The compilation target type. Can be production, profile or debug
  (defaults to production). Please note that not all targets are
  available for all systems.

For example to configure the code to configure the code to compile with
gfortran under Mac OS X in 1d do:

` $ ./configure -s macosx.gfortran -d 1`

## System configuration files

System dependent options can be configured using system configuration
files. These files are found in the config directory of the OSIRIS
source root. Users are encouraged to copy a configuration file that most
closely matches the system they will compiling for and editing it for
their specific needs.

These files are written using the [GNU Makefile
syntax](http://www.gnu.org/software/make/), and will be included by the
main Makefile prior to compilation.

### Minimal Options

The minimal options are as follows:

- **F90** - Set this parameter to the Fortran 9X compiler to use. This
  may include compiler options, or they can be set using compilation
  targets as described below. This can be set to a wrapper compiler
  (e.g. mpif90) that will automatically include the MPI and/or HDF5
  libraries.
- **C** - Set this parameter to the C compiler to use. As in the Fortran
  compiler, this may include compiler options, or they can be set using
  compilation targets as described below.
- **UNDERSCORE** - Set this to the external name convention of the
  Fortran compiler. Available options are *FORTRANNOUNDERSCORE*,
  *FORTRANSINGLEUNDERSCORE* (the most common) and
  *FORTRANDOUBLEUNDERSCORE*.
- **MPI_FCOMPILEFLAGS**, **MPI_FLINKFLAGS** - Required Fortran flags for
  compiling / linking MPI. If no flags are required (e.g. using mpif90
  for the Fortran compiler) don't set these.
- **H5_FCOMPILEFLAGS**, **H5_FLINKFLAGS** - Required Fortran flags for
  compiling / linking HDF5. If no flags are required (e.g. using h5pfc
  for the Fortran compiler) don't set these.

Here's an example of a minimal system configuration file using gfortran
for compilation:

` # Fortran compiler and flags`
` F90 = gfortran -fno-range-check -ffree-line-length-none`
` `
` # C Compiler`
` cc  = gcc`
` `
` # Fortran external name convention`
` UNDERSCORE = FORTRANSINGLEUNDERSCORE`
` `
` # MPI`
` MPI_FCOMPILEFLAGS = $(shell /opt/openmpi/1.6-gfortran/bin/mpif77 --showme:compile)`
` MPI_FLINKFLAGS    = $(shell /opt/openmpi/1.6-gfortran/bin/mpif77 --showme:link)`
` `
` # HDF5`
` H5_FCOMPILEFLAGS = -I/opt/hdf5/1.8.9-ompi-gfortran/include`
` H5_FLINKFLAGS    = -L/opt/hdf5/1.8.9-ompi-gfortran/lib -lhdf5_fortran -lhdf5 -lz -lm \`
`                    -L/opt/szip/2.1/lib  -lsz`

### Code Options

OSIRIS can be compiled using only the options specified in the previous
section, which will use the default compile options. For additional
configuration, the following parameters control compile time options:

- **PRECISION** - This controls the code precision and can be set to
  either SINGLE or DOUBLE. This overrides any settings in os-param.f90.
  If not set the code will use the definitions in os-param.f90.
- **SIMD** - This enables the use of SIMD/hardware specific code. Set
  this parameter to the appropriate architecture to use e.g. SSE for
  x86/SSE systems. If not set the code will use the default Fortran
  code. Currently available options are SSE (x86/SSE extensions), AVX
  (x86/AVX extensions), BGP (IBM BlueGene/P), and BGQ (IBM BlueGene/Q).
- **PARALLEL_IO** - Set this to __PARALLEL_IO__ to use parallel I/O
  when doing diagnostics. This option requires that HDF5 has parallel
  I/O enabled. Please note that in some systems/runs the serial
  algorithm (which has been heavily optimized) may have better
  performance.
- **TIMER** - Set this parameter to use a system specific timer. If not
  specified the code will default to using __POSIX_TIMER__ which is
  based on gettimeofday() that has a 1us resolution.

Here's an example using double precision, SSE code, not using parallel
I/O and using the Mach (OS X) timer routines:

` # Algorithm precision`
` # Set to DOUBLE or SINGLE. If not set the code will use the definition in os-param.f90`
` PRECISION = DOUBLE`
` `
` # SIMD code`
` # Set to the required SIMD algorithm (e.g. SSE or BGP). If not set the code`
` # will only compile the standard FORTRAN algorithm.`
` SIMD = SSE`
` `
` # Parallel I/O`
` # If the HDF5 library has parallel I/O support uncomment the following to use it.`
` # Note: in some systems/runs the serial I/O will be more efficient.`
` #PARALLEL_IO = __PARALLEL_IO__`
` `
` # Profiling timer`
` # Set the following to a system optimized timer, or leave commented to use a default one.`
` # This only affects code profiling measurements.`
` TIMER = __MACH_TIMER__`

### Compilation Targets

The code allows for (optionally) defining three different sets of flags
to quickly change the type of compilation targets when configuring the
code, without requiring different configuration files. These are:

- Production - This is the normal set of compilation flags to be used
  for production runs. Set the *F90FLAGS_production* parameter for this
  purpose.
- Profile - These are the compilation flags to be used for code
  profiling. Set the *F90FLAGS_profile* parameter for this purpose.
- Debug - This option will turn on debugging flags, e.g. array bounds
  checking, for development purposes. Set the *F90FLAGS_debug* parameter
  for this purpose.

Here's an example of compilation target flags for gfortran for an Intel
Core i7 CPU running Mac OS X that will use the 'Instruments'
applications for profiling:

` F90FLAGS_production = -pipe -Ofast -march=corei7 `
` F90FLAGS_profile    = -g $(F90FLAGS_production) `
` F90FLAGS_debug      = -pipe -g -fbacktrace -fbounds-check \`
`                       -Wall -fimplicit-none -pedantic \`
`                       -Wimplicit-interface -Wconversion  -Wsurprising \`
`                       -Wunderflow -ffpe-trap=invalid,zero,overflow`

### Advanced Options

Instructions for enabling OpenMP, using MPE, PAPI and SIONlib will be
published shortly.

## Setting the binary name

By default, the OSIRIS binaries will be named osiris-rXXX-YD.e, where
XXX is the SVN revision number, and Y is the number of dimensions of the
binary (1 through 3). This can be changed using the VERSION option, for
example when SVN is not installed in some system, or the user wants to
add additional information to the binary name (e.g. the mpi library or
compiler used)

Here's an example of compilation version flag that prepends 'bgq-' to
the file version:

`  VERSION = bgq-r$(shell svnversion -c ../source/ | awk -F: '{ print $$2 }' )`

This produced a file name like "osiris-bgq-r466-3D.e"

## Compiling the code

After configuring OSIRIS the code can be compiled by just calling make
(or gmake on AIX systems)

` $ make`

After compilation you will have 2 new files in the ./bin directory: i)
an OSIRIS binary with a long name describing the particular compilation
options used (e.g. *osiris-64-ompi-debug-develop-2.2.362M-2D.e*) and ii)
a symbolic link to the former with a short name just describing the
dimensions of the compiled binary (e.g. *osiris-2D.e*). The make command
can also be used to quickly change the dimensions and/or the compilation
target type:

` $ make [1d|2d|3d] [production|profile|debug]`

These options will be saved so if you do for example "make 2d", the next
time you call "make" you will be compiling for 2D.