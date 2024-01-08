# This is meant to be included by the osiris Makefile


# Add local path to vpath
vpath %.f03 ioperf

# Object list for the ioperf tool
objs_ioperf = $(subst os-main.o,ioperf.o,$(objs))

# Rule for the ioperf tool
ioperf : $(objs_ioperf)
	@echo "[LINK] $(@F)"
	@cd $(builddir) && $(F90)  $(F90F) -I$(builddir) -o $@ $(notdir $(objs_ioperf)) $(LDF)
	@chmod 755 $(builddir)/$@
	@echo "Copying binary to ioperf-$(VERSION)-$(OS_DIM)D.$(BIN_EXT)"
	@cp $(builddir)/$@ $(bindir)/ioperf-$(VERSION)-$(OS_DIM)D.$(BIN_EXT)
	@echo "Creating symbolic link"
	@ln -f -s $(bindir)/ioperf-$(VERSION)-$(OS_DIM)D.$(BIN_EXT) $(bindir)/ioperf-$(OS_DIM)D.$(BIN_EXT)
	@echo "Done!"

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
