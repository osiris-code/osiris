# This is meant to be included by the osiris Makefile

# grid object source files
src +=  os-grid-define.f03 os-grid-parallel.f90 \
        os-grid.f90

# Add local path to vpath
vpath %.f90 grid
vpath %.f03 grid


# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
