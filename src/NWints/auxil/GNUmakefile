# $Id: GNUmakefile,v 1.37 2005-02-21 05:12:57 edo Exp $

   OBJ_OPTIMIZE = int_norm.o nmcoef.o int_norm_2c.o int_chk_sum.o int_spcart.o 
            OBJ = do_chk_sum.o intb_chk_sum.o xlmcoeff.o xlm_norm.o int_spcart2.o

ifeq ($(TARGET),cray-sv2)
   OBJ += nmcoef2c.o 
else
   OBJ_OPTIMIZE += nmcoef2c.o 
endif
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

