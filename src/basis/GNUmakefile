include ../config/makefile.h

       LIBRARY = libutil.a
#      LIBRARY = libbasis.a
       HEADERS = bas.fh basP.fh bas_staticP.fh geobasmapP.fh basdeclsP.fh
           OBJ = basis.o basisP.o bas_input.o
   LIB_TARGETS = testbasis testbasis.o

ifeq ($(TARGET),DELTA)
   LIB_DEFINES = -DBASIS_LIBRARY="'library'"
else
   LIB_DEFINES = -DBASIS_LIBRARY="'$(SRCDIR)/basis/library'"
endif

     USES_BLAS = basis.F basisP.F

include ../config/makelib.h

testbasis:	testbasis.o $(LIBRARY) 
	$(FC) $(FFLAGS) -o $@ testbasis.o $(LIBS)


basis.o:	basP.fh	geobasmapP.fh basdeclsP.fh
basP.fh:	bas_staticP.fh
	@touch basP.fh
