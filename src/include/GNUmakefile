#$Id: GNUmakefile,v 1.14 1999-06-02 04:19:12 d3e129 Exp $
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
	@echo Nothing to be done

# In this directory, clean does not do anything -- we only want get rid of
# up the include files if there's a serious cleanup.

clean:
	
cleanF:
	
depend:
	
source:
	
realclean:
	$(RM) *.h *.fh
	$(RM) *~ \#*\#
