# $Id: GNUmakefile,v 1.10 2000-06-05 17:38:33 edo Exp $

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
                 hnd_initd.o  \
                 hnd_d1emem.o  \
                 hnd_stvintd.o \
                 hnd_dderij.o \
                 hnd_sd_ij.o  \
                 hnd_tvd_ij.o \
                 hnd_ddtxyz.o \
                 hnd_hlf_ij.o \
                 hnd_dvxyz.o  \
                 hnd_initdd.o \
                 hnd_dd1emem.o \
                 hnd_stvintdd.o \
                 hnd_hlfd_ij.o \
                 hnd_ddvxyz.o \
                 hnd_deriaj.o \
                 hnd_to_nw_order.o \
                 hnd_terminate.o 

       HEADERS = hnd_whermt.fh \
                 hnd_rys.fh \
                 hnd_tol.fh

       LIBRARY = libnwints.a

     USES_BLAS = hnd_stvint.F hnd_stvintd.F hnd_stvintdd.F
include ../../config/makelib.h

