# $Id: GNUmakefile,v 1.1.1.1 1994-03-29 06:44:37 d3g681 Exp $

          OBJ = inp.o
      LIBRARY = libinp.a
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o input.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

