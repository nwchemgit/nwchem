# $Id: GNUmakefile,v 1.2 1994-08-22 03:10:01 d3g681 Exp $

      LIBRARY = libutil.a
#     LIBRARY = libinp.a
          OBJ = inp.o
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o input.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

