# This is meant to be included by the osiris Makefile

# ZDF source files
src += zdf.c zdf-interface.f03 zdf-parallel.c zdf-par-interface.f03

# Add local path to vpath
vpath %.c zdf
vpath %.h zdf
vpath %.f03 zdf

# Specify dependencies
$(builddir)/zdf.o: zdf.c zdf.h

$(builddir)/zdf-interface.o: zdf-interface.f03 zdf.c zdf.h

$(builddir)/zdf-parallel.o: zdf-parallel.c zdf-parallel.h zdf.h

$(builddir)/zdf-par-interface.o: zdf-par-interface.f03 zdf-parallel.c zdf-parallel.h zdf.h


# If not at top level go there and make
ifeq ("$(toplevel)","")
all:
	$(MAKE) -C ..
clean:
	$(MAKE) -C .. clean
endif