#$Id: GNUmakefile,v 1.12 1995-02-03 18:21:22 d3h449 Exp $

include ../config/makefile.h

       LIBRARY = libblas.a

#
# Determine if single or double precision is required
#

     VPATH = double:single

ifeq ($(TARGET),SGITFP)
       PRECISION = double
endif
ifeq ($(TARGET),SGI)
       PRECISION = double
endif
ifeq ($(TARGET),SUN)
       PRECISION = double
endif
ifeq ($(TARGET),IBM)
       PRECISION = double
endif
ifeq ($(TARGET),SP1)
       PRECISION = double
endif
ifeq ($(TARGET),KSR)
       PRECISION = single   
endif
ifeq ($(TARGET),CRAY-T3D)
       PRECISION = single
endif
ifeq ($(TARGET),DELTA)
       PRECISION = double
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
OBJ_OPTIMIZE =  ddot.o    dnrm2.o    dspr.o    idamax.o \
                dasum.o   dgemm.o    drot.o    dswap.o    lsame.o \
                daxpy.o   dgemv.o    dscal.o   dtrmm.o    xerbla.o \
                dcopy.o   dger.o     dspmv.o   dtrmv.o    dtrsm.o dsyrk.o
else
OBJ_OPTIMIZE =  sdot.o    snrm2.o    sspr.o    isamax.o \
                sasum.o   sgemm.o    srot.o    sswap.o    lsame.o \
                saxpy.o   sgemv.o    sscal.o   strmm.o    xerbla.o \
                scopy.o   sger.o     sspmv.o   strmv.o    strsm.o ssyrk.o
endif
ifeq ($(TARGET),CRAY-T3D)
OBJ_OPTIMIZE = idamax.o
endif


include ../config/makelib.h


