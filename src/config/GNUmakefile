#
# $Id: GNUmakefile,v 1.15 2008-06-20 20:19:30 d3p307 Exp $
#

include makefile.h

  DEPENDCFLAGS = 

ifeq ($(TARGET),CRAY-T3D)
# Assume here that are cross compiling on a YMP
  DEPENDCFLAGS = -Tcray-ymp 
endif

ifeq ($(TARGET),cray-sv2)
# Needed to avoid having depend.x compiled as an application.
# With this flag depend.x will run on the OS node.
# From Ted Packwood
  DEPENDCFLAGS = -hcommand
endif

ifeq ($(TARGET),SP)
# Must use plain cc not mpcc
  CC = cc
endif

ifeq ($(TARGET),LAPI)
# LAPI different from SP now
# Must use plain cc not mpcc
  CC = cc
endif

ifeq ($(TARGET),LAPI64)
# LAPI different from SP now
# Must use plain cc not mpcc
  CC = cc
endif

ifeq ($(TARGET),HPUX)
  DEPENDCFLAGS = -Aa -D_HPUX_SOURCE +e
endif

ifeq ($(TARGET),LINUX)
  CC = gcc
endif

ifeq ($(TARGET),BGL)
 CC=gcc
endif

ifeq ($(TARGET),BGP)
 CC=gcc
endif

$(BINDIR)/depend.x:	depend.c
	$(CC) $(DEPENDCFLAGS) -g -o $@ $^
	/bin/rm -f depend.o

clean:
	/bin/rm -f $(BINDIR)/depend.x depend.o

realclean:	clean
	/bin/rm -f nwchem_config.h NWCHEM_CONFIG


