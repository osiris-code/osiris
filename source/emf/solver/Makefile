# This is meant to be included by the osiris Makefile

# emf object source files
src +=  os-emf-solver-yee.f03 os-emf-solver-stencil.f03 os-emf-solver-ndfx.f03 \
        os-emf-solver-lehe.f03 os-emf-solver-fei.f03

src += os-emf-solver.f03

# Add local path to vpath
vpath %.f03 emf/solver

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
