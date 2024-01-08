# This is meant to be included by the osiris Makefile

# emf object source files
src +=  os-emf-diag-define.f03 os-emf-define.f03 os-emf-gridcenter.f90

os-emf-define.f03: os-emf-solver-define.f03 os-emf-bnd-define.f03

# boundary conditions
include emf/boundary/Makefile

# field solver
include emf/solver/Makefile

# diagnostics
include emf/diagnostics/Makefile

src += os-emf-es-solver.f03 os-emf-interpolate.f03 \
  os-emf-gridval.f90 os-emf-ncifilter.f90 os-emf.f03

# Add local path to vpath
vpath %.f03 emf
vpath %.f90 emf

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
