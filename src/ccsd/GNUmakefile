#     $Id$

include ../config/makefile.h

ifdef USE_IMAX_OPENMP_TRPDRV
SUBDIRS += module
endif

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
     moints_trp.F \
     ccsd_zvecs.F \
     ccsd_z2pm.F \
     ccsd_tzinta.F \
     ccsd_tzintb.F \
     ccsd_zsig1.F \
     ccsd_zsig2.F \
     ccden_interm2.F \
     ccden_1pdm.F \
     geaxpy.F \
     ccsd_mktau_ga.F \
     ccden_2pdmb.F \
     aoccsd2.F \
     ccsd_fsig1.F \
     ccsd_fsig2.F \
     ccsd_trpdrv_bgp2.F \
     ccsd_trpdrv_offload.F \
     ccsd_trpdrv_openacc.F \
     ccsd_trpdrv_openmp_imax.F \
     moints_trp.F

ifeq ($(TARGET),BGP)
  OBJ_OPTIMIZE += ccsd_trpdrv_bgp2.o ccsd_tengy_bgp2.o ccsd_tengy_bgp.o
  LIB_DEFINES += -DBGP
endif

ifdef USE_MIC_TRPDRV
  OBJ_OPTIMIZE += ccsd_trpdrv_offload.o
  LIB_DEFINES += -DUSE_MIC_TRPDRV
endif
ifeq ($(_FC),xlf)
  LIB_DEFINES += -DTRPMIXP_OFF
endif
ifdef EXTERNAL_GA_PATH
  GA_WAPI.H=${EXTERNAL_GA_PATH}/include/ga-wapi.h
else
  GA_WAPI.H=${NWCHEM_TOP}/src/tools/build/global/src/ga-wapi.h
endif
HAVE_SET_GA_PROPERTY_DISABLES = $(shell cat  ${GA_WAPI.H} | awk ' /wnga_set_property/ {print "Y"}')
ifeq ($(HAVE_SET_GA_PROPERTY),Y)
      LIB_DEFINES += -DHAVE_SET_GA_PROPERTY
endif

ifdef USE_IMAX_OPENMP_TRPDRV

  OBJ_OPTIMIZE += ccsd_trpdrv_openmp_imax.o

  OBJ_OPTIMIZE += ccsd_trpdrv_omp_reduce_f.o

  FOPTIONS += -O3 -fiopenmp  -fopenmp-targets=spir64="-mllvm -vpo-paropt-opt-data-sharing-for-reduction=false -mllvm -vpo-paropt-atomic-free-reduction-par-global=false" -switch offload_modvars -mllvm -vpo-paropt-atomic-free-reduction-slm=true -qmkl -DMKL_ILP64 -I"${MKLROOT}/include" -I ${NWCHEM_TOP}/src/ccsd/module -mllvm -vpo-paropt-dispatch-codegen-version=1 -switch -use-host-usm-for-implicit-reduction-map

  COPTIONS:=$(filter-out -fopenmp,$(COPTIONS))
  COPTIONS:=$(filter-out -O1,$(COPTIONS))

  COPTIONS += -O3 -fiopenmp  -fopenmp-targets=spir64="-mllvm -vpo-paropt-opt-data-sharing-for-reduction=false -mllvm -vpo-paropt-atomic-free-reduction-par-global=false" -mllvm -vpo-paropt-atomic-free-reduction-slm=true -qmkl -DMKL_ILP64 -I"${MKLROOT}/include" -mllvm -vpo-paropt-dispatch-codegen-version=1

endif

ifdef USE_OPENACC_TRPDRV
  OBJ_OPTIMIZE += ccsd_trpdrv_openacc.o
  FOPTIONS += -DUSE_OPENACC_TRPDRV
  ifeq ($(_FC),pgf90)
      FOPTIONS += -Mextend -acc -cuda -cudalib=cublas
  endif
  ifeq ($(_FC),gfortran)
      FOPTIONS += -ffree-form -fopenacc -lcublas
  endif
endif


ifeq ($(ARMCI_NETWORK),MPI-PR)
   LIB_DEFINES += -DACC_STRIPS
endif

   LIB_DEFINES +=-DUSE_F90INTERFACE

 LIBRARY = libccsd.a

include ../config/makelib.h

