#$Id: GNUmakefile,v 1.16 1995-11-03 20:41:39 d3g681 Exp $

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
ifeq ($(TARGET),DECOSF)
       PRECISION = double
endif
ifeq ($(TARGET),LINUX)
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
                dcopy.o   dger.o     dspmv.o   dtrmv.o    dtrsm.o \
                dsyrk.o   dsyr2k.o   dsymm.o   dtrsv.o    dsyr2.o \
		dsymv.o 
else
OBJ_OPTIMIZE =  sdot.o    snrm2.o    sspr.o    isamax.o \
                sasum.o   sgemm.o    srot.o    sswap.o    lsame.o \
                saxpy.o   sgemv.o    sscal.o   strmm.o    xerbla.o \
                scopy.o   sger.o     sspmv.o   strmv.o    strsm.o \
                ssyrk.o   ssyr2k.o   ssymm.o   strsv.o    ssyr2.o \
		ssymv.o 
endif
ifeq ($(TARGET),CRAY-T3D)
OBJ_OPTIMIZE = idamax.o
endif


include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
