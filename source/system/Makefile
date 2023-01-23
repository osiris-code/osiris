# This is meant to be included by the osiris Makefile

# species object source files
src += os-sys-multi-c.c os-sys-multi.f03

# Add local path to vpath
vpath %.f03 system
vpath %.c system

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif