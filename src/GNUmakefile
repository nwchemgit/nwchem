
# $Id: GNUmakefile,v 1.17 1994-11-10 23:47:48 gg502 Exp $

include config/makefile.h

  SUBDIRS = $(NWSUBDIRS)

#
# This should be the first target so that a simple make builds
# the entire world
#

.PRECIOUS:	$(BINDIR)/nwchem nwchem
$(BINDIR)/nwchem nwchem:	directories libraries
	$(MAKE) nwchem.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(BINDIR)/nwchem nwchem.o $(LIBS)
	$(INSTALL)
#
# This to link only and avoid tedious subdir makes
#

link nwchem_link:    nwchem.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $(BINDIR)/nwchem nwchem.o $(LIBS)
	$(INSTALL)
#
# This for linking test version
#
nwchem_test:	
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ nwchem.o $(LIBS)
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
	-$(RM) *.db movecs 

sngl_to_dbl dbl_to_sngl:
	for dir in $(SUBDIRS); do \
		$(MAKE)	-C $$dir $@ || exit 1 ;  \
	done
        
#
# This rule to make necessary directories
#

directories:
	test -d $(LIBDIR) || $(MKDIR) $(LIBDIR)
	test -d $(BINDIR) || $(MKDIR) $(BINDIR)
	


