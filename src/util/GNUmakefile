# $Id: GNUmakefile,v 1.5 1994-04-26 19:42:06 d3g681 Exp $

include ../config/makefile.h

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
else
  EXTRA_OBJ = 
endif

       BLAS = dfill.o ifill.o rsg.o

        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o $(BLAS) $(EXTRA_OBJ) 
    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh

    USES_BLAS = util.fh

include ../config/makelib.h

