#$Id: GNUmakefile,v 1.31 2002-01-19 03:04:14 zzhang Exp $

include ../config/makefile.h

       LIBRARY = libblas.a

#
# Determine if single or double precision is required
#

     VPATH = double:single

       PRECISION = double
  
ifeq ($(TARGET),KSR)
       PRECISION = single   
endif
ifeq ($(TARGET),CRAY-T3D)
       PRECISION = single
endif
ifeq ($(TARGET),CRAY-T3E)
       PRECISION = single
endif


ifndef PRECISION
error:
	@echo BLAS makefile does not have PRECISION defined
	@exit 1
endif


#
# Use full vectorization compiler options
#

ifdef FVECTORIZE
       FOPTIMIZE := $(FVECTORIZE)
endif

ifeq ($(PRECISION),double)
#######################################################################
#               ########     ########     ########     ########     
OBJ_OPTIMIZE =  dasum.o      daxpy.o      dcopy.o      ddot.o         \
                dgemm.o      dgemv.o      dger.o       dnrm2.o        \
                drot.o       dscal.o      dspmv.o      dspr.o         \
                dswap.o      dsymm.o      dsymv.o      dsyr2.o        \
                dsyr2k.o     dsyrk.o      dtrmm.o      dtrmv.o        \
                dtrsm.o      dtrsv.o      dznrm2.o     idamax.o       \
		izamax.o \
                lsame.o      zaxpy.o      dcabs1.o     zcopy.o        \
                zdotc.o      zdscal.o     zgemm.o      zgemv.o        \
                zgerc.o      zhemm.o      zhemv.o      zher2.o        \
                zher2k.o     zherk.o      zscal.o      zswap.o        \
                ztrmm.o      ztrmv.o      ztrsm.o      ztrsv.o        \
                xerbla.o     zdotu.o   zgebak.o  zgehrd.o          \
		zgetrs.o  zladiv.o	zlanhs.o  zlarft.o  zlassq.o  \
		zung2r.o  zgebal.o  zgesv.o   zhseqr.o  zlahqr.o zlarf.o \
		zlarfx.o  zlaswp.o  zunghr.o  zgeev.o   zgetf2.o zlacgv.o \
		zlahrd.o	zlarfb.o  zlascl.o  zlatrs.o  zungqr.o \
		zgehd2.o  zgetrf.o  zlacpy.o  zlange.o	zlarfg.o  zlaset.o  \
		ztrevc.o izamax.o zgeru.o dzasum.o dlabad.o\


else
#######################################################################
#               ########     ########     ########     ########     
OBJ_OPTIMIZE =  sasum.o      saxpy.o      scopy.o      sdot.o         \
                sgemm.o      sgemv.o      sger.o       snrm2.o        \
                srot.o       sscal.o      sspmv.o      sspr.o         \
                sswap.o      ssymm.o      ssymv.o      ssyr2.o        \
                ssyr2k.o     ssyrk.o      strmm.o      strmv.o        \
                strsm.o      strsv.o      scnrm2.o     isamax.o       \
                lsame.o      caxpy.o                   ccopy.o        \
                cdotc.o      csscal.o     cgemm.o      cgemv.o        \
                cgerc.o      chemm.o      chemv.o      cher2.o        \
                cher2k.o     cherk.o      cscal.o      cswap.o        \
                ctrmm.o      ctrmv.o      ctrsm.o      ctrsv.o        \
                xerbla.o
endif
ifeq ($(TARGET),CRAY-T3D)
OBJ_OPTIMIZE = idamax.o
endif
ifeq ($(TARGET),CRAY-T3E)
OBJ_OPTIMIZE = idamax.o
endif

ifeq ($(TARGET),SP)
ifdef USE_ESSL
OBJ_OPTIMIZE = lsame.o xerbla.o
endif
endif
ifeq ($(TARGET),LAPI)
ifdef USE_ESSL
OBJ_OPTIMIZE = lsame.o xerbla.o dcabs1.o
endif
endif



include ../config/makelib.h


