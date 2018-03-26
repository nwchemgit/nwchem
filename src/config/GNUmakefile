#
# $Id$
#

include makefile.h

DEPENDCFLAGS = 
DEPEND_CC=$(CC)

ifeq ($(TARGET),HPUX)
  DEPENDCFLAGS = -Aa -D_HPUX_SOURCE +e
endif

#
# force serial compiler in case CC=mpicc used (even though it should not be)
#
ifeq ($(TARGET),LINUX)
  DEPEND_CC = gcc
endif
ifeq ($(TARGET),LINUX64)
  DEPEND_CC = gcc
endif
ifeq ($(TARGET),BGQ)
  DEPEND_CC = /usr/bin/gcc
endif


#
# use serial compilers which run on the front-end
#
ifeq ($(TARGET),LAPI)
  DEPEND_CC = cc
endif
ifeq ($(TARGET),LAPI64)
  DEPEND_CC = cc
endif
ifeq ($(TARGET),BGL)
 DEPEND_CC=gcc
endif
ifeq ($(TARGET),BGP)
 DEPEND_CC=gcc
endif

$(BINDIR)/depend.x:	depend.c
	$(DEPEND_CC) $(DEPENDCFLAGS) -g -o $@ $^
	/bin/rm -f depend.o

clean:
	/bin/rm -f $(BINDIR)/depend.x depend.o

realclean:	clean
	/bin/rm -f nwchem_config.h NWCHEM_CONFIG


