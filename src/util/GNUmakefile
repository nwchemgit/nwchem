# $Id: GNUmakefile,v 1.4 1994-04-26 04:54:52 d3e129 Exp $

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

include ../config/makelib.h

