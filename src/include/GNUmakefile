#
# This directory is a central repository for all include
# files.  The makefile in each subdirectory should contain
# a rule that keeps this directory up to date
#

include ../config/makefile.h

includes:
	for dir in $(SUBDIRS); do \
		echo Making include_stamp in $(SRCDIR)/$$dir ; \
		(cd $(SRCDIR)/$$dir; $(MAKE) include_stamp) ; \
	done

include_stamp:
	echo Nothing to be done

realclean clean:
	/bin/rm -f *.h *.fh
	/bin/rm -f *~ \#*\#
