# $Id: GNUmakefile,v 1.14 1994-08-21 22:00:27 d3e129 Exp $

include ../../config/makefile.h
API_OBJ = int_init.o exact_mem.o int_mem.o \
          int_terminate.o int_chk_sh.o int_2e4c.o \
          int_l2e4c.o int_canon.o int_chk_init.o \
          int_1eall.o int_l1eall.o int_1eh1.o \
          int_lgen1e.o int_1eov.o cando_sp.o \
          int_1eke.o int_1epe.o int_pgen1e.o \
          int_2e3c.o int_2e2c.o int_l2e3c.o \
          int_l2e2c.o int_1e3ov.o int_l1e3ov.o \
          int_hf1sp.o int_nint.o \
          intd_1eov.o intd_1eh1.o intd_2e4c.o \
          intd_init.o exactd_mem.o int_acc.o

api_obj:	$(API_OBJ)

clean:
	/bin/rm -f $(API_OBJ)

exact_mem.o:	apiP.fh
int_chk_init.o:	apiP.fh
int_init.o:	apiP.fh	apispP.fh
int_mem.o:	apiP.fh
int_terminate.o:	apiP.fh
int_hf1sp.o:	apispP.fh

# ARFLAGS = rv
# LIBRARY = ../libnwints.a
# LIBOBJ = $(patsubst %, $(LIBRARY)(%), $(API_OBJ))
#
# $(LIBRARY):	$(LIBOBJ)
