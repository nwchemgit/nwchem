# $Id: GNUmakefile,v 1.22 1995-01-10 06:06:26 d3e129 Exp $

    OBJ = int_init.o exact_mem.o int_mem.o \
          int_terminate.o int_chk_sh.o int_2e4c.o \
          int_l2e4c.o int_canon.o int_chk_init.o \
          int_1eall.o int_l1eall.o int_1eh1.o \
          int_lgen1e.o int_1eov.o cando_sp.o \
          int_1eke.o int_1epe.o int_pgen1e.o \
          int_2e3c.o int_2e2c.o int_l2e3c.o \
          int_l2e2c.o int_1e3ov.o int_l1e3ov.o \
          int_hf1sp.o int_nint.o \
          intd_1eov.o intd_1eh1.o intd_2e4c.o \
          intd_init.o exactd_mem.o int_acc.o \
          int_1cg.o intd_2e2c.o intd_2e3c.o \
          int_mpole.o int_projpole.o \
          intp_1eov.o intp_1eke.o intp_txyz.o \
          intp_1epe.o

   LIBRARY = libnwints.a

 USES_BLAS = int_hf1sp.F   int_l1e3ov.F   int_l1eall.F \
             int_l2e2c.F   int_l2e3c.F    int_l2e4c.F \
             int_lgen1e.F  int_pgen1e.F   intd_1eh1.F \
             intd_2e4c.F   intp_txyz.F    intd_2e2c.F \
             intp_1eke.F   intp_1eov.F

include ../../config/makefile.h
include ../../config/makelib.h

exact_mem.o:	apiP.fh
int_chk_init.o:	apiP.fh
int_init.o:	apiP.fh	apispP.fh
int_mem.o:	apiP.fh
int_terminate.o:	apiP.fh
int_hf1sp.o:	apispP.fh

intd_init.o int_init.o:	int_tim.fh
