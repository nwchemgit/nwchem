# $Id$

include ../../config/makefile.h

      LIB_DEFINES = -DDBLEPREC_RYS

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
                 hnd_dervxyz.o  \
                 hnd_initdd.o \
                 hnd_dd1emem.o \
                 hnd_stvintdd.o \
                 hnd_hlfd_ij.o \
                 hnd_deriaj.o \
                 hnd_terminate.o \
                 hnd_denxyz.o \
                 hnd_elfder.o \
                 hnd_giaoh01.o \
                 hnd_giaoh11.o \
                 hnd_giaol10.o \
                 hnd_giaos10.o \
                 hnd_giaotv10.o \
                 hnd_giaobq10.o \
                 hnd_giao_init.o \
                 hnd_giao_1emem.o \
                 hnd_giao_2emem.o \
                 hnd_giagnm.o \
                 hnd_giaspd.o \
                 hnd_jkbcdf.o \
                 hnd_oechrg.o \
                 hnd_xyzgia.o \
                 hnd_giahnd.o \
                 hnd_giaxyz.o \
                 hnd_jkwrys.o \
                 hnd_spdgia.o \
                 hnd_pso.o \
                 hnd_ops.o \
                 hnd_dso.o \
                 hnd_wleg.o \
                 hnd_herxyz.o \
                 hnd_init_1eelec.o \
                 hnd_init_dso.o \
                 hnd_init_1emtpl.o \
                 hnd_init_pso.o \
                 hnd_init_ops.o \
                 hnd_veloc.o \
                 hnd_angmom.o \
                 hnd_dipxyz.o \
                 hnd_init_dip.o \
                 hnd_giaos100.o

  OBJ = hnd_vstat.o \
        
       HEADERS = hnd_whermt.fh \
                 hnd_rys.fh \
                 hnd_tol.fh \
                 hnd_pointers.fh \
                 hnd_giao.fh \
                 hnd_wleg.fh

       LIBRARY = libnwints.a

     USES_BLAS = hnd_stvint.F hnd_stvintd.F hnd_stvintdd.F\
hnd_elfder.F hnd_pso.F hnd_ops.F hnd_dso.F hnd_giaoh11.F hnd_giaoh01.F hnd_giaotv10.F hnd_giaobq10.F \
hnd_giaol10.F hnd_giaos10.F hnd_giaos100.F hnd_vstat.F hnd_veloc.F hnd_dipder.F hnd_angmom.F \
hnd_dipvel.F hnd_giahnd.F

include ../../config/makelib.h
  ifeq ($(FC),efc)
            FOPTIONS += -save
#           FOPTIMIZE += -O0 -g
  endif

