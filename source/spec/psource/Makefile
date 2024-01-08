# This is meant to be included by the osiris/spec Makefile

# partticle source object source files
src += os-psource-std.f03 os-psource-constq.f03 os-psource-beam.f03 \
       os-psource-file.f03 os-psource.f03 

# Add local path to vpath
vpath %.f03 spec/psource

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
