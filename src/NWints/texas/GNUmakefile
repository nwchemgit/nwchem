# $Id: GNUmakefile,v 1.50 2004-05-14 15:55:12 edo Exp $
 
include ../../config/makefile.h

USE_TEXAS = YEP

ifeq ($(USE_TEXAS),YEP)

ifeq ($(TARGET),SOLARIS)
ifeq ($(FC), frt)
 FOPTIMIZE += -xamtfer,tfer_i
else
 FOPTIMIZE += -inline=amtfer,tfer_i
endif
endif
ifeq ($(TARGET),SP)
 FOPTIMIZE += -Q+amtfer:tfer_i
endif
ifeq ($(TARGET),LAPI)
 FOPTIMIZE += -Q+amtfer:tfer_i
endif

      LIB_DEFINES = -DUSE_TEXAS
  OBJ = abcd.o
  OBJ_OPTIMIZE = texas_face.o \
                 texas_hf.o \
                 big_service.o \
                 getmem.o \
                 amshift.o \
                 assemblx.o \
                 spec_block.o \
                 cshneg.o \
                 d_datlog.o \
                 spec_calcint.o \
                 spec_destiny.o \
                 memoha.o \
                 precalc2.o \
                 prepint.o \
                 lab_reorder.o \
                 rt123.o \
                 shells.o \
                 time.o \
                 trobsa.o \
                 zeroint.o \
		 missing.o \
		 texas_stubs.o \
                 uniq_pairs.o \
                 transfor.o \
                 gencon.o \
                 derivat.o \
                 convert.o \
                 test_calc2.o \
                 screened_op.o
else
      OBJ = texas_stubs.o 
      OBJ_OPTIMIZE = 
endif

      USES_BLAS = amshift.f \
                  big_service.F \
                  trobsa.F  getmem.F

       LIBRARY = libnwints.a
       HEADERS = screened_op.fh

include ../../config/makelib.h
ifeq ($(FC),efc)
        FDEBUG=
	FOPTIONS+= $(FOPTIMIZE) -mP3OPT_ecg_fp_control_safety=TRUE  -g -save
else
        FOPTIONS+= $(FOPTIMIZE)
endif

txs:
	@touch texas_stubs.F
	$(MAKE) "USE_TEXAS=YEP"

removeobj:
	ar d  $(LIBRARY_PATH) $(OBJ) $(OBJ_OPTIMIZE)

