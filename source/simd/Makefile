# This is meant to be included by the osiris Makefile

ifneq ("$(SIMD)","")

# SIMD source files
src += os-spec-current-c.c os-spec-push-c.c 

# Add local path to vpath
vpath %.c simd
vpath %.h simd

# Specify dependencies
simd_s := $(shell echo $(SIMD) | tr A-Z a-z)

# Except for BGQ we use the 'd' suffix for double precision code
ifneq ("$(SIMD)","BGQ")
  ifeq ("$(PRECISION)","DOUBLE")
	simd_sp := $(simd_s)d
  else
	simd_sp := $(simd_s)
  endif
else
  simd_sp := $(simd_s)
endif

$(builddir)/os-spec-push-c.o: os-spec-push-c.c os-spec-push-$(simd_sp).c \
   os-spec-push-$(simd_sp).h vector-$(simd_s).h splines-$(simd_s).h \
   os-spec-current-$(simd_sp).h split-vec.h fortran.h

$(builddir)/os-spec-current-c.o: os-spec-current-c.c os-spec-current-$(simd_sp).c \
  os-spec-current-$(simd_sp).h vector-$(simd_s).h os-spec-push-$(simd_sp).h \
  splines-$(simd_s).h fortran.h

endif

# If not at top level go there and make
ifeq ("$(toplevel)","")
all:
	$(MAKE) -C ..
clean:
	$(MAKE) -C .. clean
endif