# This is meant to be included by the osiris Makefile

src += os-diagfile.f03 

output_formats_src = os-diagfile-zdf.f03

ifdef H5_FCOMPILEFLAGS
  output_formats_src += os-dutil-hdf5.f90 os-diagfile-hdf5.f03
  FPP += -DHDF5
  ifdef H5_HAVE_PARALLEL
	FPP += -DH5_HAVE_PARALLEL
  endif
endif

output_formats_src += os-dutil.f03

# Add local path to vpath
vpath %.f90 dutil
vpath %.f03 dutil


# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
