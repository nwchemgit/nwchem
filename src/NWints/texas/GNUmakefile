# $Id: GNUmakefile,v 1.31 1997-03-26 00:49:48 d3e129 Exp $
 
USE_TEXAS = YEP

ifeq ($(USE_TEXAS),YEP)
      LIB_DEFINES = -DUSE_TEXAS
  OBJ_OPTIMIZE = texas_face.o\
                 texas_hf.o\
                 big_service.o\
                 getmem.o\
                 amshift.o\
                 assemblx.o\
                 spec_block.o\
                 cshneg.o\
                 d_datlog.o\
                 spec_calcint.o\
                 spec_destiny.o\
                 memoha.o\
                 pak.o\
                 precalc2.o\
                 prepint.o\
                 lab_reorder.o\
                 rt123.o\
                 shells.o\
                 time.o\
                 trobsa.o\
                 zeroint.o \
		 missing.o \
		 texas_stubs.o \
                 uniq_pairs.o \
                 transfor.o\
                 gencon.o 
#                test_calc2.o
else
      OBJ = texas_stubs.o
      OBJ_OPTIMIZE = 
endif

      USES_BLAS = big_service.f

       LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h

txs:
	@touch texas_stubs.F
	$(MAKE) "USE_TEXAS=YEP"

removeobj:
	ar d  $(LIBRARY_PATH) $(OBJ) $(OBJ_OPTIMIZE)

