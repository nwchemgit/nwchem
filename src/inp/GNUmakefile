# $Id: GNUmakefile,v 1.5 1995-04-17 23:39:13 gg502 Exp $

      LIBRARY = libutil.a
#     LIBRARY = libinp.a
 OBJ_OPTIMIZE = inp.o inp_irange.o inp_ilist.o
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o inp.o inp_irange.o inp_ilist.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

