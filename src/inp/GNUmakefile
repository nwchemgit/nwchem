# $Id: GNUmakefile,v 1.6 1995-07-16 02:57:12 gg502 Exp $

      LIBRARY = libutil.a
#     LIBRARY = libinp.a
 OBJ_OPTIMIZE = inp.o inp_irange.o
ifeq ($(NWCHEM_TARGET),SGITFP)
         OBJ  = inp_ilist.o
else
 OBJ_OPTIMIZE += inp_ilist.o
endif
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o inp.o inp_irange.o inp_ilist.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

