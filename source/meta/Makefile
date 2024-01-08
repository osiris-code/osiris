# This is meant to be included by the osiris Makefile

# Meta algorithm source files
src += os-emf-pgc.f03 os-spec-pgc.f03 os-neutral-pgc.f03 \
       os-particles-pgc.f03 os-simulation-pgc.f03

# Add local path to vpath
vpath %.f03 meta

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
