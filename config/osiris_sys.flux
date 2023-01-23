##########################################################################################
# System specific configuration for OSIRIS
#   System:    Linux / Hoffman and Dawson2 clusters at UCLA
#   Compilers: Intel family ( icc, ifort ) 
##########################################################################################

##########################################################################################
# OSIRIS configuration

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __POSIX_TIMER__

# SIMD
# Uncomment one of the following lines to enable SSE / AVX optimized code 
SIMD = SSE
#SIMD = AVX

SIMD_FLAGS = -msse4.2

# Numeric precision (SINGLE|DOUBLE)
#PRECISION = SINGLE
PRECISION = DOUBLE

##########################################################################################
# Compilers

F90 = ifort -msse4.2 
F03 = ifort -stand f03 -free
FPP = gcc -C -E -x assembler-with-cpp 
cc  = icc -restrict -ansi-alias

F90FLAGS_all = -fpp -diag-disable 5117 $(SIMD_FLAGS)

# Fortran preprocessor
FPP = gcc -C -E -x assembler-with-cpp

# This flag supresses some hyper-vigilant compilier warnings that have been deemed harmless. The list of warnings
# can be found is ./source/config.mk.warnings. If you are having some strange issues that you need to debug
# comment out DISABLE_PARANOIA and/or read the ./source/config.warnings file to allow the warnings if you think
# they may help you find the issue.
DISABLE_PARANOIA = YES

##########################################################################################
# Fortran flags

# External name mangling
UNDERSCORE = FORTRANSINGLEUNDERSCORE

# Flag to enable compilation of .f03 files (so far only intel compiler requires this)
F03_EXTENSION_FLAG = -free -Tf

# Enable OpenMP code
# FPP += -D_OPENMP
# F90FLAGS_all += -openmp

# ------------------------------- Compilation Targets ------------------------------------

# Debug
F90FLAGS_debug      = $(F90FLAGS_all) -O0 -debug all -traceback \
                      -warn all -fpe-all=0
# -ftrapuv 

# Production
F90FLAGS_production = $(F90FLAGS_all) -O3 -ipo -no-prec-div

# Profile Shark
# you cannot source profile with -ipo -static (which are turned on by -fast)
# After compilation you must generate the symbol table manually (the system gets confused 
# because of the extra preprocessor, the C code does not require this)
# go to the directory where the binary is located and run:
# % dsymutil <binary-name>

F90FLAGS_profile    = $(F90FLAGS_all) -O3 -debug all 

##########################################################################################
# C flags

# -wd,981,1418,2259 \
#debug
CFLAGS_debug = -O0 -g -debug full -traceback -fp-trap-all=common -Wall -w3 -wd,981,1418,2259\
               -pedantic -std=c99 

# profile
CFLAGS_profile = -O3 $(SIMD_FLAGS) -std=c99 -debug all 

# production
# CFLAGS_production = -fast
CFLAGS_production = -ipo -O3 $(SIMD_FLAGS) -no-prec-div -ansi-alias

##########################################################################################
# Linker flags

# None required
#LDFLAGS =

##########################################################################################
# Libraries

# MPI
MPI_FCOMPILEFLAGS = $(shell mpif77 --showme:compile)
MPI_FLINKFLAGS    = $(shell mpif77 --showme:link)

# MPI supports MPI_IN_PLACE reduction
FPP += -D__HAS_MPI_IN_PLACE__

# HDF5
H5_FCOMPILEFLAGS = -I$(HDF5_INCLUDE)
H5_FLINKFLAGS    = -L$(HDF5_LIB) \
                   -lhdf5_fortran -lhdf5 -lz -lm

# Uncomment the following if HDF5 supports parallel MPIIO
# H5_HAVE_PARALLEL = 1
