include makefile.h

  DEPENDCFLAGS = 

ifeq ($(TARGET),CRAY-T3D)
# Assume here that are cross compiling on a YMP
  DEPENDCFLAGS = -Tcray-ymp 
endif

depend.x:	depend.c
	$(CC) $(DEPENDCFLAGS) -g -o $@ $^

realclean clean:
	/bin/rm -f depend.x depend.o

