# $Id: GNUmakefile,v 1.4 1994-04-07 18:56:56 d3e129 Exp $

include ../../config/makefile.h
API_OBJ = int_init.o exact_mem.o int_mem.o int_terminate.o int_chk_sh.o \
          int_2e4c.o int_l2e4c.o int_canon.o int_chk_init.o int_1eall.o \
          int_l1eall.o int_1eh1.o int_lgen1e.o int_1eov.o cando_sp.o


apiobj:	$(API_OBJ)

clean:
	/bin/rm -f $(API_OBJ)

exact_mem.o:	apiP.fh
int_chk_init.o:	apiP.fh
int_init.o:	apiP.fh
int_mem.o:	apiP.fh
int_terminate.o:	apiP.fh
