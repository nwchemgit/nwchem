#
# $Id: GNUmakefile,v 1.10 2002-01-23 18:50:03 bert Exp $
#

include makefile.h

  DEPENDCFLAGS = 

ifeq ($(TARGET),CRAY-T3D)
# Assume here that are cross compiling on a YMP
  DEPENDCFLAGS = -Tcray-ymp 
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

$(BINDIR)/depend.x:	depend.c
	$(CC) $(DEPENDCFLAGS) -g -o $@ $^
	/bin/rm -f depend.o

clean:
	/bin/rm -f $(BINDIR)/depend.x depend.o

realclean:	clean
	/bin/rm -f nwchem_config.h NWCHEM_CONFIG


