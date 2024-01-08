# This is meant to be included by the osiris Makefile

# emf diagnostics object source files
src +=  os-emf-poynting.f90 os-emf-psi.f90 os-emf-diagnostics.f03

# Add local path to vpath
vpath %.f03 emf/diagnostics
vpath %.f90 emf/diagnostics

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
