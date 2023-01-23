##########################################################################################
# System specific configuration for OSIRIS
#   System:    Mac OS X / darwin
#   Compilers: GNU family ( gcc, gfortran )
##########################################################################################

##########################################################################################
# OSIRIS configuration

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __MACH_TIMER__

# SIMD
# Uncomment the following line to use SIMD optimized code. You can use SIMD = SSE, AVX or
# AVX2. When using AVX | AVX2 you must add -mavx | -mavx2 to the compiler options
#SIMD = SSE
#SIMD = AVX
SIMD = AVX2

# PRECISION
PRECISION = SINGLE
#PRECISION = DOUBLE

##########################################################################################
# Compilers

# When compiling AVX code the following must be used to
# i) enable AVX / AVX2 code generation (-mavx / -mavx2)
# and
# ii) use the clang integrated assembler instead of the GNU based system assembler (-Wa,-q).

F90 = gfortran-9.2.0 -mavx2 -Wa,-q

# Use this to specify specific options for compiling .f03 files
F03 = $(F90)
#F03 = $(F90) -std=f2003


cc  = gcc-9.2.0 -mavx2 -Wa,-q
CC  = gcc-9.2.0 -mavx2 -Wa,-q

# Fortran preprocessor
FPP = gcc-9.2.0 -C -E -x assembler-with-cpp -D__HAS_MPI_IN_PLACE__

# This flag supresses some hyper-vigilant compilier warnings that have been deemed harmless. The list of warnings
# can be found is ./source/config.mk.warnings. If you are having some strange issues that you need to debug
# comment out DISABLE_PARANOIA and/or read the ./source/config.warnings file to allow the warnings if you think
# they may help you find the issue.
DISABLE_PARANOIA = YES

##########################################################################################
# Fortran flags

# External name mangling
UNDERSCORE = FORTRANSINGLEUNDERSCORE

# Flag to enable compilation of .f03 files (not needed for gfortran)
# F03_EXTENSION_FLAG =

# ------------------------------- Compilation Targets ------------------------------------

# -fno-range-check is required because of the random module.
# gfortran has a bug that considers -2147483648 to be outside the valid
# int32 range

# -pipe makes gfortran use pipes for internal process communication (instead of files)
#       which speeds up the ocmpilation process significantly

# -ffree-line-length-none removes all constraints on line size

F90FLAGS_all = -pipe -ffree-line-length-none -fno-range-check

# OpenMP Support
F90FLAGS_all += --openmp
FPP += -D_OPENMP

# Production
# Intel Core i7 flags

F90FLAGS_production = $(F90FLAGS_all) -Ofast -march=native

# Debug

# -std=f95 is too picky

F90FLAGS_debug      = $(F90FLAGS_all) -g -Og -fbacktrace -fbounds-check \
                      -Wall -fimplicit-none -pedantic \
                      -Wimplicit-interface -Wconversion  -Wsurprising \
                      -Wunderflow  -ffpe-trap=invalid,zero,overflow

#-ffpe-trap=underflow,denormal is usually too picky. For example, it may raise an exception when
# converting double precision values to single precision for diagnostics output if the
# value is outside of valid single precision range (e.g. 1e-40). It may also raise
# an exception if a wave amplitude gets very low in a PML region. Note that in these
# situations the value is (correctly) rounded to 0.

# Profile with Shark
F90FLAGS_profile    = -g $(F90FLAGS_production)


##########################################################################################
# C flags

CFLAGS_production = -Ofast -march=native -std=c99

CFLAGS_debug      = -Og -g -Wall -pedantic -march=native -std=c99
#-fsanitize=address

CFLAGS_profile    = -g $(CFLAGS_production)

##########################################################################################
# Linker flags

# None required
#LDFLAGS =

##########################################################################################
# Libraries

# MPI
OPENMPI_ROOT = /opt/openmpi/4.0.2

MPI_CCOMPILEFLAGS = $(shell $(OPENMPI_ROOT)/bin/mpicc --showme:compile)
MPI_FCOMPILEFLAGS = $(shell $(OPENMPI_ROOT)/bin/mpif90 --showme:compile)
MPI_FLINKFLAGS    = $(shell $(OPENMPI_ROOT)/bin/mpif90 --showme:link)

# HDF5
H5_ROOT = /opt/hdf5/1.10.5

H5_FCOMPILEFLAGS = -I$(H5_ROOT)/include
H5_FLINKFLAGS    = -L$(H5_ROOT)/lib -lhdf5_fortran -lhdf5 -lz -lm

# Uncomment the following if HDF5 supports parallel MPIIO
# H5_HAVE_PARALLEL = 1
