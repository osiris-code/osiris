# This is meant to be included by the osiris Makefile

src += os-spec-rawdiag.f03 os-spec-phasespace.f03 os-spec-diagnostics.f03

# Add local path to vpath
vpath %.f03 spec/diagnostics

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif