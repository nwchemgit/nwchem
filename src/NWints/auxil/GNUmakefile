# $Id: GNUmakefile,v 1.31 1997-02-17 20:51:17 d3g681 Exp $

   OBJ_OPTIMIZE = int_norm.o nmcoef.o nmcoef2c.o int_norm_2c.o int_chk_sum.o 
            OBJ = do_chk_sum.o intb_chk_sum.o xlmcoeff.o int_spcart.o xlm_norm.o
LIBRARY = libnwints.a
LIB_TARGETS = spherical spherical.o testnat testnat.o

USES_BLAS = int_spcart.F spherical.F xlm_norm.F xlmcoeff.F spcartP.fh

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
testnat:	testnat.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testnat.o

