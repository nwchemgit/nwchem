#     $Id$

include ../../config/makefile.h

ifdef USE_IMAX_OPENMP_TRPDRV
 OBJ_OPTIMIZE = \

  OBJ_OPTIMIZE  += ccsd_trpdrv_mkl_module.o


  FOPTIONS += -fiopenmp -fopenmp-targets=spir64="-mllvm -vpo-paropt-enable-64bit-opencl-atomics=true -mllvm -vpo-paropt-opt-data-sharing-for-reduction=false" -qmkl -DMKL_ILP64 -I"${MKLROOT}/include" -fpp -fixed -free

#OBJ_OPTIMIZE += ccsd_trpdrv_mkl_module.o
#.SECONDARY: ccsd_trpdrv_mkl_module.o
#ccsd_trpdrv_mkl_module.o: ccsd_trpdrv_mkl_module.F
#	$(FC) $(MK_MODULE_FLAG) -c $< -o $@

#-include $(ccsd_trpdrv_mkl_module_)
#.PRECIOUS: ccsd_trpdrv_mkl_module.F
#	$(FC) $(MK_MODULE_FLAG) -c $< -o $@
#USES_BLAS += ccsd_trpdrv_mkl_module.F
#ccsd_trpdrv_mkl_module.o : ccsd_trpdrv_mkl_module.F
#	$(FC) $(MK_MODULE_FLAG) -c $< -o $@

endif

  LIB_DEFINES +=-DUSE_F90INTERFACE

  LIBRARY = libccsd.a

include ../../config/makelib.h

