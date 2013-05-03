# $Id$
 
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
                 zeroint.o \
                 missing.o \
                 texas_stubs.o \
                 uniq_pairs.o \
                 transfor.o \
                 gencon.o \
                 derivat.o \
                 convert.o \
                 test_calc2.o \
                 screened_op.o \
                 trobsa.o
else
      OBJ = texas_stubs.o 
      OBJ_OPTIMIZE = 
endif

      USES_BLAS = amshift.F \
                  big_service.F \
                  trobsa.F  getmem.F

       LIBRARY = libnwints.a
       HEADERS = screened_op.fh 

ifdef XLF11
  DEFINES += -DXLF11
endif

include ../../config/makelib.h
ifeq ($(_CPU),ia64)
ifeq ($(FC),ifort)
        FDEBUG=
	FOPTIONS+= $(FOPTIMIZE) -mP3OPT_ecg_fp_control_safety=TRUE  -g -save
endif
endif

txs:
	@touch texas_stubs.F
	$(MAKE) "USE_TEXAS=YEP"

removeobj:
	ar d  $(LIBRARY_PATH) $(OBJ) $(OBJ_OPTIMIZE)

