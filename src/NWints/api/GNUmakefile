# $Id: GNUmakefile,v 1.73 1999-08-20 18:45:46 d3e129 Exp $


ifeq ($(USE_TEXAS),YEP)
      LIB_DEFINES = -DUSE_TEXAS
endif

    OBJ = cando_sp.o     exact_mem.o    exactd_mem.o \
          int_1cg.o      int_1e3ov.o    int_1eall.o \
          int_1eh1.o     int_1eke.o     int_1eov.o \
          int_1epe.o     int_2e2c.o     int_2e3c.o \
          int_acc.o      int_canon.o \
          int_chk_init.o int_chk_sh.o   int_func_c.o\
          int_init.o     int_l1e3ov.o   int_l1eall.o \
          int_l1eh1.o    int_l1eke.o    int_l1eov.o \
          int_l1epe.o    int_l2e2c.o    int_l2e3c.o \
          int_lgen1e.o   int_mem.o      int_mpolel.o \
          int_mpole.o    int_nint.o     int_pgen1e.o \
          int_projpole.o int_term.o \
          intd_1eh1.o    intd_1eov.o    intd_1eke.o \
          int_1er.o \
          intd_2e2c.o \
          intd_2e3c.o    intd_2e4c.o    intd_init.o \
          intd_1e3ov.o   intbd_init4c.o intbd_2e4c.o \
          intd_mpolel.o \
          intdd_1eh1.o   intdd_1eov.o   intdd_2e2c.o \
          intdd_2e4c.o   intdd_2e3c.o   intdd_init.o \
          intp_1eke.o    intp_1eov.o    intp_1epe.o \
          intp_2e2c.o    intp_2e3c.o    intp_mpole.o \
          intp_1e1cpe.o  intp_2e3ct.o   intp_mpolel.o \
          intp_txyz.o \
          intpd_1e1cpe.o intpd_1eke.o   intpd_1eov.o \
          intpd_1epe.o   intpd_2e2c.o   intpd_2e3c.o \
          intpd_2e3ct.o  intpd_mpolel.o \
	  intso_1e.o     

    OBJ_OPTIMIZE = int_hf1sp.o int_l2e4c.o intb_init4c.o intb_2e4c.o \
                   int_2e4c.o

    HEADERS = int_tim.fh apiP.fh int_nbf.fh
   LIBRARY = libnwints.a

 USES_BLAS = int_1e3ov.F int_hf1sp.F int_lgen1e.F intb_2e4c.F intbd_2e4c.F intd_1e3ov.F \
             intd_1eh1.F intd_2e2c.F intd_2e4c.F intd_mpolel.F intp_txyz.F intpd_mpolel.F


include ../../config/makefile.h
include ../../config/makelib.h

justo:	$(OBJ)	$(OBJ_OPTIMIZE)
	@echo "just objects made"
noo:
	rm -f $(OBJ) $(OBJ_OPTIMIZE)

txs:
	@touch cando_sp.F int_2e2c.F int_2e3c.F
	$(MAKE) "USE_TEXAS=YEP"

doc:
	rm -f integral_api.tex
	$(BINDIR)/seetex integral_api.th integral_api.tex
	rm -f int_api_spec.tex
	$(BINDIR)/seetex int_api_spec.th int_api_spec_init.th        int_api_spec.tex
	$(BINDIR)/seetex int_init.F intd_init.F int_term.F int_acc.F int_api_spec.tex
	$(BINDIR)/seetex int_api_spec_mem.th                         int_api_spec.tex
	$(BINDIR)/seetex int_mem.F                                   int_api_spec.tex
	$(BINDIR)/seetex int_api_spec_ints.th                        int_api_spec.tex
	$(BINDIR)/seetex int_1cg.F   int_1e3ov.F int_l1e3ov.F int_1eall.F int_l1eall.F  int_pgen1e.F int_api_spec.tex
	$(BINDIR)/seetex int_1eh1.F  int_l1eh1.F int_1eke.F   int_l1eke.F int_1eov.F    int_l1eov.F  int_api_spec.tex
	$(BINDIR)/seetex int_1epe.F  int_l1epe.F int_lgen1e.F int_2e2c.F  int_l2e2c.F                int_api_spec.tex
	$(BINDIR)/seetex int_l2e3c.F int_2e3c.F  int_2e4c.F   int_l2e4c.F intb_init4c.F intb_2e4c.F  int_api_spec.tex
	$(BINDIR)/seetex int_api_spec_prop.th                        int_api_spec.tex
	$(BINDIR)/seetex int_mpole.F int_projpole.F                  int_api_spec.tex
	$(BINDIR)/seetex int_api_spec_misc.th                        int_api_spec.tex
	$(BINDIR)/seetex exact_mem.F exactd_mem.F int_canon.F int_chk_init.F int_chk_sh.F            int_api_spec.tex
	$(BINDIR)/seetex int_func_c.F int_hf1sp.F int_nint.F                                         int_api_spec.tex
	cp integral_api.tex int_api_spec.tex ${TOPDIR}/doc/prog

