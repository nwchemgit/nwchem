#     $Id$

include ../config/makefile.h

 #OBJ = \

 OBJ_OPTIMIZE = \
        ccden_driver.o \
        ccden_interm2.o \
        ccden_1pdm.o \
        ccsd_driver.o \
        ccsd_sortmo.o \
        geaxpy.o \
        ccden_2pdma.o \
        ccden_2pdmb.o \
        ccsd_mktau_ga.o \
        ccden_interm4.o \
        aoccsd2.o \
        ccsd_itm.o \
        ccsd_t2pm.o \
        ccsd_pampt2.o \
        ccsd_pampt3.o \
        ccsd_trpdrv.o \
        ccsd_trpdrv_nb.o \
        ccsd_trpdrv_bgp.o \
        moints_trp.o \
        ccsd.o \
        ccsd_aux.o \
        ccsd_fsig1.o \
        ccsd_fsig2.o \
        ccsd_get_param.o \
        ccsd_iterdrv2.o \
        ccsd_mktask.o \
        ccsd_mkt2.o \
        ccsd_pdiis.o \
        ccsd_trpmo.o \
        ccsd_rdtrpo.o \
        qstat.o \
        cputime.o \
        rddcom.o \
        ccsd_input.o \
        ccsd_graddrv.o \
        ccsd_zvecs.o \
        ccsd_pzamp.o \
        ccsd_zsig1.o \
        ccsd_zsig2.o \
        ccsd_z2pm.o \
        ccsd_tzinta.o \
        ccsd_tzintb.o \
        ccsd_mem.o \
        ccsd_pstat.o \
        f_memset.o

ifdef XLF10
  OBJ += ccsd_aosig.o
else
  OBJ_OPTIMIZE += ccsd_aosig.o
endif

 USES_BLAS = \
     ccsd_itm.F \
     ccsd_t2pm.F \
     ccsd_fsig1.F \
     ccsd_fsig2.F \
     ccsd_iterdrv2.F \
     ccsd_pampt2.F \
     ccsd_pampt3.F \
     ccsd_pdiis.F \
     ccsd_sortmo.F \
     ccsd_trpdrv.F \
     ccsd_trpdrv_nb.F \
     ccsd_trpdrv_bgp.F \
     moints_trp.F \
     ccsd_zvecs.F \
     ccsd_z2pm.F \
     ccsd_tzinta.F \
     ccsd_tzintb.F \
     ccsd_zsig1.F \
     ccsd_zsig2.F \
     ccden_interm2.F \
     ccden_1pdm.F \
     geaxpy.f \
     ccsd_mktau_ga.F \
     ccden_2pdmb.F \
     aoccsd2.F \
     moints_trp.F

ifeq ($(TARGET),BGP)
  OBJ_OPTIMIZE += ccsd_trpdrv_bgp2.o ccsd_tengy_bgp2.o
  USES_BLAS += ccsd_trpdrv_bgp2.F
  LIB_DEFINES += -DBGP
endif

 #LIB_DEFINES = 
 #LIB_DEFINES = -DDEBUG_PRINT

 LIBRARY = libccsd.a

include ../config/makelib.h

 #FOPTIMIZE=
