SHELL = /bin/sh

# configfile is set in /source/Makefile to the path of the configurations file
#   (usually ../config/osiris_config, but different if doing and out-of-source build)
include $(configfile)

# default preprocessor
ifeq ("$(FPP)","")
  FPP = gcc -C -E -x assembler-with-cpp
endif

# Check Fortran compiler definitions
ifeq ("$(F90)","")
    $(error Configure script did not set F90 variable (F9x compiler). Please correct your configuration file)
endif

ifeq ("$(F03)","")
    F03 = $(F90)
endif

# Check compiler definitions
ifeq ("$(cc)","")
    $(warning Configure script did not set cc variable (C compiler), defaulting to gcc)
    cc = gcc
endif


# Get git version ID
ifeq ("$(VERSION)","")
VERSION := $(shell git --no-pager describe --tags --always --dirty)
endif


ifeq ("$(BRANCH_NAME)","")
BRANCH_NAME := $(shell git rev-parse --abbrev-ref HEAD)
endif

# Append branch name to executable and symlink if desired
ifdef APPEND_BRANCH
  ifeq ("$(BRANCH_SUFFIX)","")
    BRANCH_SUFFIX := -$(BRANCH_NAME)
  endif
endif


# Get build date
ifeq ("$(DATE)","")
DATE := $(shell date)
endif


# C Compiler flags
CF   = $(CFLAGS_$(COMPILATION_TYPE)) $(MPI_CCOMPILEFLAGS)

# Add Fortran underscore definition
CF += -D$(UNDERSCORE)

F90F = $(F90FLAGS_$(COMPILATION_TYPE)) $(MPI_FCOMPILEFLAGS) $(H5_FCOMPILEFLAGS)

FPPF = -DP_X_DIM=$(OS_DIM) -DOS_REV=\"$(VERSION)\" -D$(UNDERSCORE)
LDF  = $(LDFLAGS) $(H5_FLINKFLAGS) $(MPI_FLINKFLAGS) $(SION_FLINKFLAGS)


# Precision flags
ifneq ("$(PRECISION)","")
  ifneq ("$(PRECISION)","SINGLE")
   ifneq ("$(PRECISION)","DOUBLE")
     $(error When set, PRECISION must be either SINGLE or DOUBLE. Please correct your configuration file)
   endif
  endif
  CF   +=-DPRECISION_$(PRECISION)
  FPPF +=-DPRECISION_$(PRECISION)
endif


# SIMD flags
ifneq ("$(SIMD)","")
  ifeq ("$(PRECISION)","")
    $(error SIMD code selected, but PRECISION was not set. Please correct your configuration file)
  endif
  CF   +=-DSIMD -DSIMD_$(SIMD)
  FPPF +=-DSIMD -DSIMD_$(SIMD)
endif

# Allocatable array extensions
ifneq ("$(ALLOC_EXT)","")
  FPPF+=-D$(ALLOC_EXT)
endif

# System timers
ifneq ("$(TIMER)","")
  CF+=-D$(TIMER)
endif

# Disable file flushing
ifneq ("$(NO_FILE_FLUSH)","")
  FPPF+=-D$(NO_FILE_FLUSH)
endif

# Enable log files
ifneq ("$(USE_LOG)","")
  FPPF+=-D$(USE_LOG)
endif

# Restart io backend
ifneq ("$(RST_IO)","")
  FPPF+=-D$(RST_IO)
endif

ifneq ("$(FFTW_ROOT)","")
  FPPF += -DFFTW_ENABLED
  F90F += $(FFTW_FCOMPILEFLAGS)
  LDF  += $(FFTW_FLINKFLAGS)
endif

# binary extension
ifeq ("$(BIN_EXT)","")
  BIN_EXT = e
endif
