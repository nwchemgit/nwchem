include makefile.h

  DEPENDCFLAGS = 

ifeq ($(TARGET),CRAY-T3D)
# Assume here that are cross compiling on a YMP
  DEPENDCFLAGS = -Tcray-ymp 
endif

ifeq ($(TARGET),SP1)
# Must use plain cc not mpcc
  CC = cc
endif

$(BINDIR)/depend.x:	depend.c
	$(CC) $(DEPENDCFLAGS) -g -o $@ $^
	/bin/rm -f depend.o

realclean clean:
	/bin/rm -f $(BINDIR)/depend.x depend.o

