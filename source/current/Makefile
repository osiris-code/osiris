# This is meant to be included by the osiris Makefile

# Electric current source files
src += os-current-diagnostics.f03 os-current-boundary.f03 os-current-define.f03

# Add local path to vpath
vpath %.f03  current

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
