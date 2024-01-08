# This is meant to be included by the osiris Makefile

# particle object source files
src +=   memory.f03

# Add local path to vpath
vpath %.f03 memory

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
