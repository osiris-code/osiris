# This is meant to be included by the osiris Makefile

# species object source files
src += os-spec-define.f03 os-spec-memory.f90 os-spec-v3dep.f90 os-spec-v3int.f90 \
       os-spec-udist.f03 os-spec-sdist.f03 os-spec-tracks.f03 os-spec-tag.f03 os-spec-comm.f03 \
       os-spec-current.f03 os-spec-charge.f90 os-spec-initfields.f90 \
       os-spec-piston.f90 os-spec-boundary.f90

# diagnostics
include spec/diagnostics/Makefile

# particle source object source files
include spec/psource/Makefile

src += os-spec-loadbalance.f90 os-spec-push.f90 os-spec-extpush.f03 os-spec-vpush.f90 os-spec-accel.f90 \
       os-spec-radcool.f90 os-species.f03 os-spec-collisions.f03


# Add local path to vpath
vpath %.f03 spec
vpath %.f90 spec

# Add extra dependencies
$(builddir)/os-spec-define.o: os-psource-def.f03 os-spec-diag-def.f03

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
