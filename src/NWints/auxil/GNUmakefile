# $Id: GNUmakefile,v 1.14 1995-10-30 20:55:37 d3e129 Exp $

   OBJ_OPTIMIZE = int_normalize.o nmcoef.o nmcoef2c.o int_norm_2c.o int_chk_sum.o \
                  xlmcoeff.o
            OBJ = do_int_chk_sum.o intb_chk_sum.o
LIBRARY = libnwints.a
LIB_TARGETS = spherical spherical.o testnat testnat.o

USES_BLAS = do_int_chk_sum.F

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
testnat:	testnat.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testnat.o

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
do_int_chk_sum.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
int_chk_sum.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
int_norm_2c.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
int_normalize.o :	 $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
int_spcart.o :	 $(INCDIR)/stdio.fh spcartP.fh $(INCDIR)/mafdecls.fh 
