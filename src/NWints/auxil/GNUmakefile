# $Id: GNUmakefile,v 1.34 1997-11-07 01:46:08 d3e129 Exp $

   OBJ_OPTIMIZE = int_norm.o nmcoef.o nmcoef2c.o int_norm_2c.o int_chk_sum.o int_spcart.o
            OBJ = do_chk_sum.o intb_chk_sum.o xlmcoeff.o xlm_norm.o 

 HEADERS = spcartP.fh


 LIBRARY = libnwints.a

 LIB_TARGETS = spherical spherical.o testnat testnat.o

 USES_BLAS = int_spcart.F spherical.F xlm_norm.F xlmcoeff.F 

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
testnat:	testnat.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testnat.o

