# This is meant to be included by the osiris Makefile

ifneq ("$(FFTW_INSTALL_DIR)","")
  $(error ERROR Specifying FFTW_INSTALL_DIR in configuration file is now deprecated. \
		Instead, specify FFTW_ROOT, FFTW_FCOMPILEFLAGS, and FFTW_FLINKFLAGS. \
		See config/osiris_sys.linux.gnu for an example)
endif

# grid object source files
# add no files for now in this inital integration stub...
# src +=  fftlib.f90 fftmod.f90
src += os-fft.f03 os-fft-new.f03

# Add local path to vpath
vpath %.f03 fft

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
