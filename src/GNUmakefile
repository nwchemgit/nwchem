
# $Id: GNUmakefile,v 1.2 1994-04-01 22:36:11 d3g681 Exp $

include config/makefile.h

nwchem:	libraries
	$(MAKE) nwchem.o
	$(FC) $(FFLAGS) -o $@ nwchem.o $(LIBS)
	$(INSTALL)

# This dependency so that includes which are made by libraries rule
# are updated before nwchem is linked

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
	


