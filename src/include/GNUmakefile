#$Id$
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

include_stamp sngl_to_dbl dbl_to_sngl 64_to_32 32_to_64:
	@echo Nothing to be done

# In this directory, clean does not do anything -- we only want get rid of
# up the include files if there's a serious cleanup.

clean:
	
cleanF:
	
cleanDEP:
	
depend:
	
source:
	
realclean:
	$(RM) *.h *.fh
	$(RM) *~ \#*\#
