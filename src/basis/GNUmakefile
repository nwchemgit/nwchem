#$Id: GNUmakefile,v 1.11 1995-02-02 23:09:31 d3g681 Exp $
include ../config/makefile.h

       LIBRARY = libutil.a
#      LIBRARY = libbasis.a
       HEADERS = bas.fh basP.fh bas_staticP.fh geobasmapP.fh basdeclsP.fh
           OBJ = basis.o basisP.o bas_input.o
ifeq ($(TARGET),CRAY-T3D)
   OBJ_OPTIMIZE = basis.o 
   OBJ =  basisP.o bas_input.o
endif
   LIB_TARGETS = testbasis testbasis.o

ifeq ($(TARGET),DELTA)
   LIB_DEFINES = -DBASIS_LIBRARY="'library'"
else
   LIB_DEFINES = -DBASIS_LIBRARY="'$(SRCDIR)/basis/library'"
endif

     USES_BLAS = basis.F basisP.F

include ../config/makelib.h

testbasis:	testbasis.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testbasis.o $(LIBS)


basis.o:	basP.fh	geobasmapP.fh basdeclsP.fh
basP.fh:	bas_staticP.fh
	@touch basP.fh
