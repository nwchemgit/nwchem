# $Id: GNUmakefile,v 1.4 1995-04-17 21:32:51 gg502 Exp $

      LIBRARY = libutil.a
#     LIBRARY = libinp.a
 OBJ_OPTIMIZE = inp.o inp_irange.o
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o inp.o inp_irange.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

