# $Id: GNUmakefile,v 1.3 1994-04-22 15:16:48 d3g681 Exp $

include ../config/makefile.h

ifeq ($(TARGET),SUN)
  EXTRA_OBJ = ieeetrap.o
else
  EXTRA_OBJ = 
endif

       BLAS = dfill.o ifill.o ddot.o dcopy.o rsg.o idamax.o dgemm.o dscal.o \
              xerbla.o lsame.o

        OBJ = output.o errquit.o ffflush.o \
              print_center.o util_flush.o $(BLAS) $(EXTRA_OBJ) 
    LIBRARY = libutil.a
    HEADERS = util.fh itri.fh

include ../config/makelib.h

