# $Id: GNUmakefile,v 1.19 1996-02-27 00:31:17 d3e129 Exp $
 
ifeq ($(USE_TEXAS),YEP)
      LIB_DEFINES = -DUSE_TEXAS
  OBJ_OPTIMIZE = texas_face.o\
                 texas_hf.o\
                 big_service.o\
                 getmem.o\
                 amshift.o\
                 assemblx.o\
                 spec_block.o\
                 spec_calcint.o\
                 cshneg.o\
                 d_datlog.o\
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
                 transfor.o
else
      OBJ = texas_stubs.o
      OBJ_OPTIMIZE = 
endif

       LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h

txs:
	@touch texas_stubs.F
	$(MAKE) "USE_TEXAS=YEP"

removeobj:
	ar d  $(LIBRARY_PATH) $(OBJ) $(OBJ_OPTIMIZE)

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(precalc2.o):	 $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(spec_block.o):	 $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(texas_face.o):	 $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/bas_ibs_sfn.fh $(INCDIR)/bas_ibs_dec.fh $(INCDIR)/bas_exndcf_sfn.fh $(INCDIR)/bas_exndcf_dec.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh $(INCDIR)/apiP.fh 
$(LIBRARY_PATH)(texas_hf.o):	 $(INCDIR)/mafdecls.fh 
