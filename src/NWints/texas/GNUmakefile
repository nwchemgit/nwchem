# $Id: GNUmakefile,v 1.37 1998-06-09 23:36:57 d3e129 Exp $
 
include ../../config/makefile.h

USE_TEXAS = YEP

ifeq ($(USE_TEXAS),YEP)

ifeq ($(TARGET),SOLARIS)
 FOPTIMIZE += -inline=amtfer,tfer_i
endif
ifeq ($(TARGET),SP1)
 FOPTIMIZE += -Q+amtfer:tfer_i
endif

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
                 gencon.o\
                 derivat.o\
                 conv24r.o \
                 test_calc2.o \
                 conbuf2_der.o destdul.o
else
      OBJ = texas_stubs.o 
      OBJ_OPTIMIZE = 
endif

      USES_BLAS = amshift.f \
                  big_service.f \
                  trobsa.f 

       LIBRARY = libnwints.a

include ../../config/makelib.h

txs:
	@touch texas_stubs.F
	$(MAKE) "USE_TEXAS=YEP"

removeobj:
	ar d  $(LIBRARY_PATH) $(OBJ) $(OBJ_OPTIMIZE)

