# This is meant to be included by the osiris Makefile

# restart object source files
src += os-restart.f90

# Add local path to vpath
vpath %.f90 restart

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

