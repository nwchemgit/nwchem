# $Id$

   OBJ_OPTIMIZE = int_norm.o nmcoef.o int_norm_2c.o int_chk_sum.o int_spcart.o 
            OBJ = do_chk_sum.o intb_chk_sum.o xlmcoeff.o xlm_norm.o int_spcart2.o

   OBJ_OPTIMIZE += nmcoef2c.o 
 HEADERS = spcartP.fh


 LIBRARY = libnwints.a

 LIB_TARGETS = spherical spherical.o testnat testnat.o

 USES_BLAS = int_spcart.F spherical.F xlm_norm.F xlmcoeff.F  int_norm_2c.F int_spcart2.F

include ../../config/makefile.h
include ../../config/makelib.h

spherical:	spherical.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ spherical.o $(LIBS)
testnat:	testnat.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testnat.o

