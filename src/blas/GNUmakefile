#$Id$
ifeq ($(BLAS_SIZE),4)
USE_I4FLAGS=y
endif
include ../config/makefile.h

       LIBRARY = libnwcblas.a

#
# USE both double and single now (short term solution for added single precision
# data in GA). Still keep the original layout, in case later need to be changed
# again.
	VPATH = double single

#       PRECISION = double single
#  
#ifeq ($(TARGET),KSR)
#       PRECISION = single   
#endif
#
#
#ifndef PRECISION
#error:
#	@echo BLAS makefile does not have PRECISION defined
#	@exit 1
#endif


#
# Use full vectorization compiler options
#

ifdef FVECTORIZE
       FOPTIMIZE := $(FVECTORIZE)
endif

# xerbla and lsame are taken out from single
#ifeq ($(PRECISION),double)
#######################################################################
#               ########     ########     ########     ########     
OBJ_OPTIMIZE =  dasum.o      daxpy.o      dcopy.o      ddot.o         \
                dgemm.o      dgemv.o      dger.o       dnrm2.o        \
                drot.o       dscal.o      dspmv.o      dspr.o         \
                dswap.o      dsymm.o      dsymv.o      dsyr2.o        \
                dsyr2k.o     dsyrk.o      dtrmm.o      dtrmv.o        \
                dtrsm.o      dtrsv.o      dznrm2.o     idamax.o       \
		izamax.o     drotg.o \
                lsame.o      zaxpy.o      dcabs1.o     zcopy.o        \
                zdrot.o      zrotg.o \
                zdotc.o      zdscal.o     zgemm.o      zgemv.o        \
                zgerc.o      zhemm.o      zhemv.o      zher2.o        \
                zher2k.o     zherk.o      zscal.o      zswap.o        \
                ztrmm.o      ztrmv.o      ztrsm.o      ztrsv.o        \
                zsyrk.o      zsyr2k.o     zsymm.o      zher.o \
                xerbla.o     zgeru.o dzasum.o dsyr.o   zdotu.o
#else
#######################################################################
#               ########     ########     ########     ########     
OBJ_OPTIMIZE +=sasum.o      saxpy.o      scopy.o      sdot.o         \
                sgemm.o      sgemv.o      sger.o       snrm2.o        \
                srot.o       sscal.o      sspmv.o      sspr.o         \
                sswap.o      ssymm.o      ssymv.o      ssyr2.o        \
                ssyr2k.o     ssyrk.o      strmm.o      strmv.o        \
                strsm.o      strsv.o      scnrm2.o     isamax.o       \
                scabs1.o \
                caxpy.o      ccopy.o      csrot.o      icamax.o  \
                cdotc.o      csscal.o     cgemm.o      cgemv.o        \
                cgerc.o      chemm.o      chemv.o      cher2.o        \
                cher2k.o     cherk.o      cscal.o      cswap.o        \
                ctrmm.o      ctrmv.o      ctrsm.o      ctrsv.o        \
                cgeru.o scasum.o cdotu.o
#endif

ifeq ($(TARGET),MACX)
  ifdef USE_VECLIB
    OBJ_OPTIMIZE = lsame.o 
  endif
endif
ifeq ($(TARGET),IBM)
  ifdef USE_ESSL
    OBJ_OPTIMIZE = lsame.o xerbla.o  izamax.o dcabs1.o
  endif
endif
ifeq ($(TARGET),IBM64)
  ifdef USE_INTEGER4
    ifdef USE_ESSL
      OBJ_OPTIMIZE = lsame.o xerbla.o  izamax.o dcabs1.o
    endif
  endif
endif
ifeq ($(TARGET),SP)
  ifdef USE_ESSL
    OBJ_OPTIMIZE = lsame.o xerbla.o izamax.o dcabs1.o
  endif
endif
ifeq ($(TARGET),LAPI)
ifdef USE_ESSL
OBJ_OPTIMIZE = lsame.o xerbla.o dcabs1.o izamax.o
endif
endif

#caching for dgemm
ifeq ($(TARGET),LINUX64)
     ifeq ($(_CPU),x86_64)
           LIB_DEFINES += -DCACHE1M
     endif
     ifeq ($(_CPU),ia64)
           LIB_DEFINES += -DCACHE6M
     endif
endif

include ../config/makelib.h


