#$Id: GNUmakefile,v 1.13 1995-11-16 20:38:09 d3g681 Exp $
include ../config/makefile.h

       LIBRARY = libutil.a
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

basP.fh:	bas_staticP.fh
	@touch basP.fh

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
bas_input.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh bas.fh 
basis.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh geobasmapP.fh $(INCDIR)/context.fh basdeclsP.fh basP.fh basP.fh 
basisP.o :	 $(INCDIR)/inp.fh geobasmapP.fh basdeclsP.fh basP.fh 
testbasis.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh bas.fh 
