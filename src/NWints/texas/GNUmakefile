# $Id: GNUmakefile,v 1.36 1998-06-05 16:04:50 d3j191 Exp $
 
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

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(texas_face.o):	$(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/bas_ibs_sfn.fh $(INCDIR)/bas_ibs_dec.fh $(INCDIR)/bas_exndcf_sfn.fh $(INCDIR)/bas_exndcf_dec.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh $(INCDIR)/apiP.fh 
$(LIBRARY_PATH)(precalc2.o):	$(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(texas_hf.o):	$(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(spec_block.o):	$(INCDIR)/mafdecls.fh 
