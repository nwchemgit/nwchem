       HEADERS = bas.fh basP.fh bas_staticP.fh geobasmapP.fh basdeclsP.fh
           OBJ = basis.o newbasis.o bas_input.o
       LIBRARY = libbasis.a
   LIB_TARGETS = testbasis testbasis.o

include ../config/makefile.h
include ../config/makelib.h

testbasis:	testbasis.o $(LIBRARY) 
	$(FC) $(FFLAGS) -o $@ testbasis.o $(LIBS)


basis.o:	basP.fh	geobasmapP.fh basdeclsP.fh
basP.fh:	bas_staticP.fh
	@touch basP.fh
