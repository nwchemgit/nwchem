# $Id: GNUmakefile,v 1.3 1994-08-27 13:36:03 d3g681 Exp $

      LIBRARY = libutil.a
#     LIBRARY = libinp.a
 OBJ_OPTIMIZE = inp.o
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o input.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

