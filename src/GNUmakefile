
# $Id: GNUmakefile,v 1.7 1994-08-22 00:01:08 d3g681 Exp $

include config/makefile.h

  SUBDIRS = $(NWSUBDIRS)

#
# This should be the first target so that a simple make builds
# the entire world
#

nwchem:	libraries
	$(MAKE) nwchem.o
	$(FC) $(FLDOPT) -o $@ nwchem.o $(LIBS)
	$(INSTALL)
#
# This to link only and avoid tedious subdir makes
#

nwchem_link:    nwchem.o
	$(FC) $(FLDOPT) -o nwchem nwchem.o $(LIBS)
	$(INSTALL)
#
# This for linking test version
#
nwchem_test:	
	$(FC) $(FLDOPT) -o nwchem_test nwchem.o $(LIBS)
	$(INSTALL)

# This dependency so that includes which are made by libraries rule
# are updated before nwchem is compiled

libraries:
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir || exit 1 ;  \
	done

clean:
	/bin/rm -f nwchem core *.o *.a *.trace
	for dir in $(SUBDIRS) ../tests; do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir clean || exit 1 ;  \
	done

realclean:	clean
	/bin/rm -f *~ \#*\#
	for dir in $(SUBDIRS) ../tests; do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir realclean || exit 1 ;  \
	done
	/bin/rm -f $(LIBDIR)/*.a $(BINDIR)/*
	/bin/rm -f *.db movecs 

sngl_to_dbl dbl_to_sngl:
	for dir in $(SUBDIRS); do \
		$(MAKE)	 $(MAKEOVERRIDES) -C $$dir $@ || exit 1 ;  \
	done
        
#
# This rule to make necessary directories
#

directories:
	if [! test -d $(LIBDIR) ]; then
		mkdir $(LIBDIR) 
	fi
	if [! test -d $(BINDIR) ]; then
		mkdir $(BINDIR) 
	fi
	


