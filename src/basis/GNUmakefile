#$Id: GNUmakefile,v 1.21 1996-05-14 17:54:31 d3g681 Exp $
include ../config/makefile.h

       LIBRARY = libutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh basdeclsP.fh geobasmapP.fh  \
                 bas_exndcf.fh bas_exndcf_dec.fh bas_exndcf_sfn.fh \
                 bas_ibs.fh bas_ibs_dec.fh bas_ibs_sfn.fh 
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

     TESTLIBS = -lutil

include ../config/makelib.h

testbasis:	testbasis.o
	$(FC) $(FFLAGS) -o $@ testbasis.o  $(LDFLAGS) $(LIBS) $(TESTLIBS)

basP.fh:	bas_staticP.fh
	@touch basP.fh

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(bas_input.o):	$(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh ./bas.fh 
$(LIBRARY_PATH)(basis.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh ./geobasmapP.fh $(INCDIR)/context.fh ./basdeclsP.fh ./bas_ibs_sfn.fh ./bas_ibs_dec.fh ./bas_exndcf_sfn.fh ./bas_exndcf_dec.fh ./bas_exndcf.fh ./basP.fh 
$(LIBRARY_PATH)(basisP.o):	$(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh ./geobasmapP.fh ./basdeclsP.fh ./bas_exndcf_sfn.fh ./bas_exndcf_dec.fh ./basP.fh 
$(LIBRARY_PATH)(bas_contrib.o):	$(INCDIR)/geom.fh ./bas.fh 
$(LIBRARY_PATH)(testbasis.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geom.fh ./geobasmapP.fh ./basdeclsP.fh ./bas_ibs_sfn.fh ./bas_ibs_dec.fh ./basP.fh ./bas.fh 
