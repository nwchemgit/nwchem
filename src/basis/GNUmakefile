#$Id: GNUmakefile,v 1.16 1995-11-28 22:24:42 og845 Exp $
include ../config/makefile.h

       LIBRARY = libutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh geobasmapP.fh basdeclsP.fh
           OBJ = bas_input.o bas_contrib.o
  OBJ_OPTIMIZE = basis.o basisP.o

ifeq ($(TARGET),CRAY-T3D)
   OBJ_OPTIMIZE = basis.o 
   OBJ =  basisP.o bas_input.o bas_contrib.o
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
$(LIBRARY_PATH)(bas_contrib.o):	 $(INCDIR)/geom.fh ./bas.fh 
$(LIBRARY_PATH)(bas_input.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh ./bas.fh 
$(LIBRARY_PATH)(basis.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh ./geobasmapP.fh $(INCDIR)/context.fh ./basdeclsP.fh ./basP.fh ./basP.fh 
$(LIBRARY_PATH)(basisP.o):	 $(INCDIR)/inp.fh ./geobasmapP.fh ./basdeclsP.fh ./basP.fh 
$(LIBRARY_PATH)(testbasis.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh ./bas.fh 
