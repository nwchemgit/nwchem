
# $Id: GNUmakefile,v 1.24 1995-06-29 16:21:39 d3g681 Exp $

include config/makefile.h

  SUBDIRS =  $(NWSUBDIRS)

  LIB_DEFINES =  -DCOMPILATION_DATE="'`date`'" \
                 -DCOMPILATION_DIR="'$(NWCHEM_TOP)'"

#
# This should be the first target so that a simple make builds
# the entire world
#

.PHONY:	nwchem_config

.PRECIOUS:	$(BINDIR)/nwchem nwchem
$(BINDIR)/nwchem nwchem:	directories libraries
	$(MAKE) nwchem.o stubs.o
	$(LINK.f) -o $(BINDIR)/nwchem nwchem.o stubs.o $(LIBS)
	/bin/rm -f nwchem.o
	$(INSTALL)
#
# This to link only and avoid tedious subdir makes
#

link nwchem_link:    nwchem.o stubs.o
	$(LINK.f) -o $(BINDIR)/nwchem nwchem.o stubs.o $(LIBS)
	/bin/rm -f nwchem.o
	$(INSTALL)
#
# This for linking test version
#
test nwchem_test:	nwchem.o stubs.o
	$(LINK.f) -o nwchem_test nwchem.o stubs.o $(LIBS)
	/bin/rm -f nwchem.o
	$(INSTALL)

# This dependency so that includes which are made by libraries rule
# are updated before nwchem is compiled

libraries:
	for dir in $(SUBDIRS); do \
		$(MAKE)	-C $$dir || exit 1 ;  \
	done

clean:
	-$(RM) nwchem core *.o *.a *.trace
	for dir in $(SUBDIRS) ../tests; do \
		$(MAKE)	-C $$dir clean || exit 1 ;  \
	done

realclean:
	-$(RM) nwchem core *.o *.a *.trace
	-$(RM) *~ \#*\#
	for dir in $(SUBDIRS) ../tests; do \
		$(MAKE)	-C $$dir realclean || exit 1 ;  \
	done
	-$(RM) $(LIBDIR)/*.a $(BINDIR)/*
	-$(RM) *.db movecs fort.41 fort.44

sngl_to_dbl dbl_to_sngl:
	for dir in $(SUBDIRS); do \
		$(MAKE)	-C $$dir $@ || exit 1 ;  \
	done
        
#
# This rule to make necessary subdirectories in lib and bin (normally)
#

directories:
	test -d $(LIBDIR) || $(MKDIR) -p $(LIBDIR)
	test -d $(BINDIR) || $(MKDIR) -p $(BINDIR)

# This target will change config/nwchem_config.h include file which
# controls what we're going to build.
nwchem_config:
ifeq ($(NWCHEM_MODULES),)
	@echo
	@echo -----
	@echo You must set NWCHEM_MODULES to the list of modules to be \
		included in NWChem builds.
	@$(MAKE) -C config -f make_nwchem_config list-options
else
	$(MAKE) -C config -f make_nwchem_config NWCHEM_MODULES="$(NWCHEM_MODULES)"
endif

# This is a convenience target that will make sure all subdirectories
# listed in $(SUBDIRS) exist. Useful because CVS will only update
# directories that already exist in your working tree, so this is a
# way ot make sure all that are relevant to the code exist.

make_source_dirs:
ifneq ($(SUBDIRS),)
	for dir in $(SUBDIRS); do \
	   test -d $$dir  || $(MKDIR) $$dir; \
	done
endif

# This is a convenience target to remove source directories that are
# configured out of the program.  Useful to save disk space.

remove_unused_dirs:
ifneq ($(EXCLUDED_SUBDIRS),)
	for dir in $(EXCLUDED_SUBDIRS); do \
	   $(RM) -r $$dir; \
	done
endif
