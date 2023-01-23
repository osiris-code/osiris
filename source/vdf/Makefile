# This is meant to be included by the osiris Makefile

# vdf and wall object source files
src += os-vdf-define.f03 os-vdf-memory.f90 os-vdf-reportfile.f03 os-vdf-comm.f90 \
       os-vdf-math.f90 os-vdf-smooth.f90 \
       os-vdf-average.f03 os-vdf.f03 os-vdf-report.f03 os-vdf-interpolate.f90 \
       os-wall-define.f90 os-wall-comm.f90 os-wall.f90 \
       os-cmplx-vdf.f03 os-cmplx-vdf-math.f03

# Add local path to vpath
vpath %.f90 vdf
vpath %.f03 vdf


# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
