# $Id: GNUmakefile,v 1.6 1999-01-26 01:06:36 d3e129 Exp $

include ../../config/makefile.h

  OBJ_OPTIMIZE = hnd_whermt.o \
                 hnd_rt123.o  \
                 hnd_root4.o  \
                 hnd_root5.o  \
                 hnd_droot.o  \
                 hnd_1emem.o  \
                 hnd_init.o   \
                 hnd_set_acc.o \
                 hnd_stvint.o  \
                 hnd_sxyz.o   \
                 hnd_txyz.o   \
                 hnd_to_nw_order.o \
                 hnd_terminate.o 

       HEADERS = hnd_whermt.fh \
                 hnd_rys.fh \
                 hnd_tol.fh

       LIBRARY = libnwints.a

     USES_BLAS = hnd_stvint.F
include ../../config/makelib.h

