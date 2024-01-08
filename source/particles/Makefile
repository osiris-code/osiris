# This is meant to be included by the osiris Makefile

# particle object source files
src +=   os-particles-define.f03 os-particles-charge.f90 os-particles.f03


# Add local path to vpath
vpath %.f03 particles
vpath %.f90 particles

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
