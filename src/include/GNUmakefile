#
# This directory is a central repository for all include
# files.  The makefile in each subdirectory should contain
# a rule that keeps this directory up to date
#

include ../config/makefile.h

SUBDIRS = $(NWSUBDIRS)

includes:
	for dir in $(SUBDIRS); do \
		echo Making include_stamp in $(SRCDIR)/$$dir ; \
		(cd $(SRCDIR)/$$dir; $(MAKE) include_stamp) ; \
	done

include_stamp sngl_to_dbl dbl_to_sngl:
	echo Nothing to be done

realclean clean:
	mv msgids.fh msgids.FH
	/bin/rm -f *.h *.fh
	/bin/rm -f *~ \#*\#
	mv msgids.FH msgids.fh
