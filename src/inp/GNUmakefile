# $Id$

      LIBRARY = libnwcutil.a
 OBJ_OPTIMIZE = inp.o inp_irange.o
          OBJ = c_inp.o hnd_rdfree.o
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
	$(LINK.f) $(FFLAGS) -o $@ $^ $(LIBS)

hnd_rdfree:	hnd_rdfree.o
	$(LINK.f) -o $@ $^ $(LIBS) $(LIBS)

