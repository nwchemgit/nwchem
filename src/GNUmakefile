
# $Id: GNUmakefile,v 1.1.1.1 1994-03-29 06:44:24 d3g681 Exp $

include config/makefile.h

nwchem:	nwchem.o libraries
	$(FC) $(FFLAGS) -o $@ nwchem.o $(LIBS)
	$(INSTALL)

# This dependency so that includes which are made by libraries rule
# are updated before nwchem is linked

nwchem:	libraries

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
	


