# This is meant to be included by the osiris Makefile

# Virtual base class for PRNG modules
src += os-random-class.f03

# Multiple random number generator modules
src += os-random-mt.f03 os-random-r250.f03 os-random-cmwc.f03 \
       os-random-mwc.f03 os-random-kiss.f03 os-random-hash.f03 \
       os-random-class-grid.f03

# General interface for creating PRNG objects  
src += os-random.f03

# Add local path to vpath
vpath %.f03 random

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif
