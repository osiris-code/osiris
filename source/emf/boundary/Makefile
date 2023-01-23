# This is meant to be included by the osiris Makefile

# emf boundary source files
src +=  os-emf-lindman.f90 os-emf-vpml.f03 \
        os-emf-pec.f90 os-emf-pmc.f90 os-emf-boundary.f90

# Add local path to vpath
vpath %.f03 emf/boundary
vpath %.f90 emf/boundary

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
