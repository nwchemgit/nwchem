
# $Id: GNUmakefile,v 1.4 1994-04-28 21:59:16 d3e129 Exp $

include config/makefile.h

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

# This dependency so that includes which are made by libraries rule
# are updated before nwchem is compiled

libraries:
	for dir in $(SUBDIRS); do \
		echo Making library in $$dir ; \
		(cd $$dir; $(MAKE)) || exit 1 ; \
	done

clean:
	/bin/rm -f nwchem core *.o *.a *.trace
	for dir in $(SUBDIRS) ../tests; do \
		echo Making clean in $$dir ; \
		(cd $$dir; $(MAKE) clean) || exit 1 ; \
	done

realclean:	clean
	/bin/rm -f *~ \#*\#
	for dir in $(SUBDIRS) ../tests; do \
		echo Making realclean in $$dir ; \
		(cd $$dir; $(MAKE) realclean) || exit 1; \
	done
	/bin/rm -f $(LIBDIR)/*.a $(BINDIR)/*
	/bin/rm -f *.db

sngl_to_dbl dbl_to_sngl:
	for dir in $(SUBDIRS); do \
		echo Making $@ in $$dir ; \
		(cd $$dir; $(MAKE) $@) || exit 1 ; \
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
	


