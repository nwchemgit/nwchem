# $Id: GNUmakefile,v 1.12 1995-05-18 07:18:32 d3e129 Exp $

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
