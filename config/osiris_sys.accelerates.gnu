##########################################################################################
# System specific configuration for OSIRIS
#   System:    Linux / Accelerates cluster
#   Compilers: GNU family ( gcc, gfortran ) 
##########################################################################################

##########################################################################################
# OSIRIS configuration

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __POSIX_TIMER__

# SIMD
# Uncomment the following line to use SIMD optimized code. You can use SIMD = SSE or
# SIMD = AVX. When using AVX you must add -mavx to the compiler options
# SIMD = SSE
SIMD = AVX

# PRECISION
PRECISION = SINGLE
#PRECISION = DOUBLE

##########################################################################################
# Compilers

# When compiling AVX code the following must be used to i) enable AVX code generation (-mavx)
F90 = gfortran -mavx
F03 = gfortran -mavx -std=gnu

cc  = gcc -Wa,-q

# Fortran preprocessor
FPP = gcc -C -nostdinc -E -x assembler-with-cpp -D__HAS_MPI_IN_PLACE__

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
F90FLAGS_debug = $(F90FLAGS_all) -g -Og -fbacktrace -fbounds-check \
                 -Wall -fimplicit-none -pedantic \
                 -Wimplicit-interface -Wconversion  -Wsurprising \
                 -Wunderflow  -ffpe-trap=invalid,zero,overflow
                      
#-ffpe-trap=underflow,denormal is usually too picky. For example, it may raise an exception when
# converting double precision values to single precision for diagnostics output if the
# value is outside of valid single precision range (e.g. 1e-40). It may also raise 
# an exception if a wave amplitude gets very low in a PML region. Note that in these 
# situations the value is (correctly) rounded to 0.

# Profile
F90FLAGS_profile = -g $(F90FLAGS_production)
##########################################################################################
# C flags
CFLAGS_production = -Ofast -march=native -std=c99
CFLAGS_debug      = -Og -g -Wall -pedantic -march=native -std=c99 -fsanitize=address
CFLAGS_profile    = -g $(CFLAGS_production) 

##########################################################################################
# Linker flags

# None required
#LDFLAGS =

##########################################################################################
# Libraries
# MPI_FCOMPILEFLAGS = -I$(MPI_INCLUDE)
# MPI_FLINKFLAGS    = -Wl,-rpath -Wl,$(MPI_LIB) -L$(MPI_LIB) \
#                     -lmpifort -Wl,--enable-new-dtags
# MPI_CCOMPILEFLAGS = -I$(MPI_INCLUDE)
# 
# H5_FCOMPILEFLAGS = -I$(HDF5_INCLUDE) 
# H5_FLINKFLAGS    = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5 -lz -ldl -lm \
#                    -Wl,-rpath -Wl,$(HDF5_LIB) -Wl,--enable-new-dtags
