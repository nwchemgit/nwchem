# $Id: GNUmakefile,v 1.16 1995-11-28 20:43:29 gg502 Exp $

      LIBRARY = libutil.a
 OBJ_OPTIMIZE = inp.o inp_irange.o
          OBJ = c_inp.o
ifeq ($(TARGET),SGITFP)
         OBJ  += inp_ilist.o
else
 OBJ_OPTIMIZE += inp_ilist.o
endif
  LIB_TARGETS = test.o test
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o inp.o inp_irange.o inp_ilist.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(c_inp.o):	 $(INCDIR)/types.f2c.h 
$(LIBRARY_PATH)(inp.o):	 ./inpP.fh ./inp.fh 
$(LIBRARY_PATH)(inp_ilist.o):	 ./inpP.fh 
$(LIBRARY_PATH)(inp_irange.o):	 ./inpP.fh 
$(LIBRARY_PATH)(test.o):	 ./inp.fh 
$(LIBRARY_PATH)(test_ilist.o):	 ./inp.fh 
