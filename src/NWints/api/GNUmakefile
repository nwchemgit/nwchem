# $Id: GNUmakefile,v 1.33 1995-12-22 01:39:12 d3e129 Exp $


ifeq ($(USE_TEXAS),YEP)
      LIB_DEFINES = -DUSE_TEXAS
endif

    OBJ = cando_sp.o     exact_mem.o    exactd_mem.o \
          int_1cg.o      int_1e3ov.o    int_1eall.o \
          int_1eh1.o     int_1eke.o     int_1eov.o \
          int_1epe.o     int_2e2c.o     int_2e3c.o \
          int_2e4c.o     int_acc.o      int_canon.o \
          int_chk_init.o int_chk_sh.o   \
          int_init.o     int_l1e3ov.o   int_l1eall.o \
          int_l1eh1.o    int_l1eke.o    int_l1eov.o \
          int_l1epe.o    int_l2e2c.o    int_l2e3c.o \
          int_l2e4c.o    int_lgen1e.o   int_mem.o \
          int_mpole.o    int_nint.o     int_pgen1e.o \
          int_projpole.o int_term.o \
          intd_1eh1.o    intd_1eov.o    intd_2e2c.o \
          intd_2e3c.o    intd_2e4c.o    intd_init.o \
          intp_1eke.o    intp_1eov.o    intp_1epe.o \
          intp_2e2c.o    intp_2e3c.o    intp_mpole.o \
          intp_1e1cpe.o  intp_txyz.o \
          intb_init4c.o  intb_2e4c.o
    OBJ_OPTIMIZE = int_hf1sp.o 

   LIBRARY = libnwints.a

 USES_BLAS = intd_1eh1.F   intd_2e2c.F    intb_2e4c.F \
             int_hf1sp.F   int_lgen1e.F   intd_2e2c.F \
             intp_txyz.F

include ../../config/makefile.h
include ../../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(cando_sp.o):	 $(INCDIR)/geobasmapP.fh ./candoP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(exact_mem.o):	 $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(exactd_mem.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(int_1cg.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh ./apiP.fh 
$(LIBRARY_PATH)(int_1e3ov.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_1eall.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_1eh1.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_1eke.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_1eov.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_1epe.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_2e2c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(int_2e3c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(int_2e4c.o):	 $(INCDIR)/pstat.fh ./int_tim.fh $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_acc.o):	 ./apiP.fh 
$(LIBRARY_PATH)(int_chk_init.o):	 ./apiP.fh 
$(LIBRARY_PATH)(int_chk_sh.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_hf1sp.o):	 $(INCDIR)/stdio.fh 
$(LIBRARY_PATH)(int_init.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh ./int_tim.fh $(INCDIR)/global.fh ./candoP.fh $(INCDIR)/bas.fh ./apispP.fh ./apiP.fh 
$(LIBRARY_PATH)(int_l1e3ov.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_l1eall.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_l1eh1.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_l1eke.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_l1eov.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_l1epe.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_l2e2c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_l2e3c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_l2e4c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_lgen1e.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_mem.o):	 $(INCDIR)/stdio.fh ./apiP.fh 
$(LIBRARY_PATH)(int_mpole.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh ./apiP.fh 
$(LIBRARY_PATH)(int_nint.o):	 $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_pgen1e.o):	 $(INCDIR)/geobasmapP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(int_projpole.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh ./apiP.fh 
$(LIBRARY_PATH)(int_term.o):	 $(INCDIR)/mafdecls.fh ./apispP.fh ./apiP.fh 
$(LIBRARY_PATH)(intb_2e4c.o):	 $(INCDIR)/stdio.fh ./numb_qP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(intb_init4c.o):	 ./numb_qP.fh 
$(LIBRARY_PATH)(intd_1eh1.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(intd_1eov.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(intd_2e2c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(intd_2e3c.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(intd_2e4c.o):	 $(INCDIR)/pstat.fh ./int_tim.fh $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(intd_init.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh ./apispP.fh ./apiP.fh 
$(LIBRARY_PATH)(intp_1e1cpe.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh ./apiP.fh 
$(LIBRARY_PATH)(intp_1eke.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(intp_1eov.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(intp_1epe.o):	 $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh 
$(LIBRARY_PATH)(intp_2e2c.o):	 $(INCDIR)/stdio.fh $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(intp_2e3c.o):	 $(INCDIR)/stdio.fh $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh $(INCDIR)/bas.fh ./apiP.fh 
$(LIBRARY_PATH)(intp_mpole.o):	 $(INCDIR)/stdio.fh $(INCDIR)/geomP.fh $(INCDIR)/geobasmapP.fh $(INCDIR)/basdeclsP.fh $(INCDIR)/basP.fh ./apiP.fh 
$(LIBRARY_PATH)(intp_txyz.o):	 $(INCDIR)/geomP.fh 
