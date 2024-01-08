# This is meant to be included by the osiris Makefile

src += os-zpulse-std.f03 os-zpulse-wall.f03 os-zpulse-point.f03 \
       os-zpulse-speckle.f03 os-zpulse.f03

# Add local path to vpath
vpath %.f03 zpulse

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
