# $Id: GNUmakefile,v 1.32 1997-02-20 10:51:06 d3e129 Exp $

   OBJ_OPTIMIZE = int_norm.o nmcoef.o nmcoef2c.o int_norm_2c.o int_chk_sum.o 
            OBJ = do_chk_sum.o intb_chk_sum.o xlmcoeff.o xlm_norm.o int_spcart.o


 LIBRARY = libnwints.a

 LIB_TARGETS = spherical spherical.o testnat testnat.o

 USES_BLAS = int_spcart.F spherical.F xlm_norm.F xlmcoeff.F spcartP.fh

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
testnat:	testnat.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testnat.o

