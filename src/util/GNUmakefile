# $Id: GNUmakefile,v 1.2 1994-04-08 20:41:07 d3g681 Exp $

include ../config/makefile.h

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
else
  EXTRA_OBJ = 
endif

        OBJ = dfill.o ifill.o output.o rsg.o errquit.o ffflush.o \
              ddot.o dcopy.o $(EXTRA_OBJ) 
    LIBRARY = libutil.a
    HEADERS = util.fh

include ../config/makelib.h

