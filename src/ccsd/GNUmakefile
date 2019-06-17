#     $Id$

include ../config/makefile.h

 OBJ_OPTIMIZE = \
        ccsd_driver.o \
        ccsd_sortmo.o \
        ccsd_mktau_ga.o \
        aoccsd2.o \
        ccsd_itm.o \
        ccsd_itm_omp.o \
        ccsd_t2pm.o \
        ccsd_pampt2.o \
        ccsd_pampt3.o \
	ccsd_idx1.o \
	ccsd_idx2.o \
	ccsd_idx34.o \
	ccsd_tengy.o \
        ccsd_trpdrv.o \
        ccsd_trpdrv_nb.o \
        ccsd_trpdrv_omp.o \
	ccsd_trpdrv_omp_mixprec.o \
        moints_trp.o \
        ccsd.o \
        ccsd_aux.o \
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
        ccsd_mem.o \
        ccsd_pstat.o \
        convert_single_double.o

ifdef XLF10
  OBJ += ccsd_aosig.o
else
  OBJ_OPTIMIZE += ccsd_aosig.o
endif

 USES_BLAS = \
     ccsd_itm.F \
     ccsd_itm_omp.F \
     ccsd_t2pm.F \
     ccsd_iterdrv2.F \
     ccsd_pampt2.F \
     ccsd_pampt3.F \
     ccsd_idx1.F \
     ccsd_idx2.F \
     ccsd_idx34.F \
     ccsd_pdiis.F \
     ccsd_sortmo.F \
     ccsd_trpdrv.F \
     ccsd_trpdrv_nb.F \
     ccsd_trpdrv_omp.F \
     ccsd_trpdrv_omp_mixprec.F \
     ccsd_mktau_ga.F \
     aoccsd2.F \
     moints_trp.F

ifeq ($(TARGET),BGP)
  OBJ_OPTIMIZE += ccsd_trpdrv_bgp2.o ccsd_tengy_bgp2.o ccsd_tengy_bgp.o
  USES_BLAS += ccsd_trpdrv_bgp2.F
  LIB_DEFINES += -DBGP
endif

ifdef USE_MIC_TRPDRV
  OBJ_OPTIMIZE += ccsd_trpdrv_offload.o
  USES_BLAS    += ccsd_trpdrv_offload.F
  LIB_DEFINES += -DUSE_MIC_TRPDRV
endif
ifeq ($(_FC),xlf)
  LIB_DEFINES += -DTRPMIXP_OFF
endif
HAVE_SET_GA_PROPERTY_DISABLES = $(shell cat ${NWCHEM_TOP}/src/tools/build/global/src/ga-wapi.h | awk ' /wnga_set_property/ {print "Y"}')
ifeq ($(HAVE_SET_GA_PROPERTY),Y)
      LIB_DEFINES += -DHAVE_SET_GA_PROPERTY
endif


ifeq ($(ARMCI_NETWORK),MPI-PR)
   LIB_DEFINES += -DACC_STRIPS
endif

   LIB_DEFINES +=-DUSE_F90INTERFACE

 LIBRARY = libccsd.a

include ../config/makelib.h

