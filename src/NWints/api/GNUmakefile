# $Id: GNUmakefile,v 1.25 1995-05-09 18:01:03 d3e129 Exp $

    OBJ = cando_sp.o     exact_mem.o    exactd_mem.o \
          int_1cg.o      int_1e3ov.o    int_1eall.o \
          int_1eh1.o     int_1eke.o     int_1eov.o \
          int_1epe.o     int_2e2c.o     int_2e3c.o \
          int_2e4c.o     int_acc.o      int_canon.o \
          int_chk_init.o int_chk_sh.o   int_hf1sp.o \
          int_init.o     int_l1e3ov.o   int_l1eall.o \
          int_l1eh1.o    int_l1eke.o    int_l1eov.o \
          int_l1epe.o    int_l2e2c.o    int_l2e3c.o \
          int_l2e4c.o    int_lgen1e.o   int_mem.o \
          int_mpole.o    int_nint.o     int_pgen1e.o \
          int_projpole.o int_terminate.o \
          intd_1eh1.o    intd_1eov.o    intd_2e2c.o \
          intd_2e3c.o    intd_2e4c.o    intd_init.o \
          intp_1eke.o    intp_1eov.o    intp_1epe.o \
          intp_2e2c.o    intp_2e3c.o    intp_mpole.o \
          intp_txyz.o \
          intb_init4c.o  intb_2e4c.o

   LIBRARY = libnwints.a

 USES_BLAS = int_hf1sp.F   int_l1e3ov.F   int_l1eall.F \
             int_l2e2c.F   int_l2e3c.F    int_l2e4c.F \
             int_lgen1e.F  int_pgen1e.F   intd_1eh1.F \
             intd_2e4c.F   intp_txyz.F    intd_2e2c.F \
             intp_1eke.F   intp_1eov.F

include ../../config/makefile.h
include ../../config/makelib.h

exact_mem.o:	apiP.fh
exactd_mem.o:	apiP.fh
int_1cg.o:	apiP.fh
int_2e2c.o:	apiP.fh
int_2e3c.o:	apiP.fh
int_2e4c.o:	int_tim.fh
int_acc.o:	apiP.fh
int_chk_init.o:	apiP.fh
int_hf1sp.o:	apispP.fh
int_init.o:	apiP.fh	int_tim.fh	apispP.fh
int_mem.o:	apiP.fh
int_mpole.o:	apiP.fh
int_projpole.o:	apiP.fh
int_terminate.o:	apiP.fh	apispP.fh
intb_2e4c.o:	numb_qP.fh
intb_init4c.o:	numb_qP.fh
intd_2e2c.o:	apiP.fh
intd_2e3c.o:	apiP.fh
intd_2e4c.o:	int_tim.fh
intd_init.o:	apiP.fh	apispP.fh
intp_2e2c.o:	apiP.fh
intp_2e3c.o:	apiP.fh
intp_mpole.o:	apiP.fh

intd_init.o int_init.o:	int_tim.fh

