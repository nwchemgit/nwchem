# $Id: GNUmakefile,v 1.20 1997-01-29 19:36:40 d3g681 Exp $

      LIBRARY = libutil.a
 OBJ_OPTIMIZE = inp.o inp_irange.o
          OBJ = c_inp.o 
#hnd_rdfree.o
ifeq ($(TARGET),SGITFP)
         OBJ  += inp_ilist.o
else
 OBJ_OPTIMIZE += inp_ilist.o
endif
  LIB_TARGETS = test.o test hnd_rdfree.o hnd_rdfree
      HEADERS = inp.fh


include ../config/makefile.h
include ../config/makelib.h

test:	test.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

hnd_rdfree:	hnd_rdfree.o
	$(LINK.f) -o $@ $^ $(LIBS) $(LIBS)

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(hnd_rdfree.o):	./inp.fh 
$(LIBRARY_PATH)(inp.o):	./inpP.fh ./inp.fh 
$(LIBRARY_PATH)(test_ilist.o):	./inp.fh 
$(LIBRARY_PATH)(test.o):	./inp.fh 
$(LIBRARY_PATH)(c_inp.o):	$(INCDIR)/types.f2c.h 
$(LIBRARY_PATH)(inp_irange.o):	./inpP.fh 
$(LIBRARY_PATH)(inp_ilist.o):	./inpP.fh 
