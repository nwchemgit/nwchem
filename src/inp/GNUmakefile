# $Id: GNUmakefile,v 1.8 1995-11-07 09:19:35 apsmith Exp $

      LIBRARY = libutil.a
#     LIBRARY = libinp.a
# OBJ_OPTIMIZE = inp.o inp_irange.o
#ifeq ($(NWCHEM_TARGET),SGITFP)
#         OBJ  = inp_ilist.o
#else
# OBJ_OPTIMIZE += inp_ilist.o
#endif
     OBJ = inp.o inp_irange.o inp_ilist.o c_inp.o
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o inp.o inp_irange.o inp_ilist.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
inp.o :	 inpP.fh inp.fh 
inp_ilist.o :	 inpP.fh 
inp_irange.o :	 inpP.fh 
test.o :	 inp.fh 
test_ilist.o :	 inp.fh 
