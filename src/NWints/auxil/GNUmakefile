# $Id: GNUmakefile,v 1.13 1995-10-17 08:52:32 d3g681 Exp $

   OBJ_OPTIMIZE = int_normalize.o nmcoef.o nmcoef2c.o int_norm_2c.o int_chk_sum.o \
                  xlmcoeff.o
            OBJ = do_int_chk_sum.o 
LIBRARY = libnwints.a
LIB_TARGETS = spherical spherical.o

USES_BLAS = do_int_chk_sum.F

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
do_int_chk_sum.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
int_chk_sum.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
int_norm_2c.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
int_normalize.o :	 $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
int_spcart.o :	 $(INCDIR)/stdio.fh spcartP.fh $(INCDIR)/mafdecls.fh 
