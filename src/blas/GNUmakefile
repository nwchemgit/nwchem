
include ../config/makefile.h

       LIBRARY = libblas.a

#
# Determine if single or double precision is required
#

       VPATH = double

ifeq ($(TARGET),KSR)
       VPATH = single   
endif
ifeq ($(TARGET),CRAY)
       VPATH = single
endif

ifeq ($(VPATH),double)
         OBJ =  ddot.o    dnrm2.o    dspr.o    idamax.o \
                dasum.o   dgemm.o    drot.o    dswap.o    lsame.o \
                daxpy.o   dgemv.o    dscal.o   dtrmm.o    xerbla.o \
                dcopy.o   dger.o     dspmv.o   dtrmv.o
else
         OBJ =  sdot.o    snrm2.o    sspr.o    isamax.o \
                sasum.o   sgemm.o    srot.o    sswap.o    lsame.o \
                saxpy.o   sgemv.o    sscal.o   strmm.o    xerbla.o \
                scopy.o   sger.o     sspmv.o   strmv.o
endif


include ../config/makelib.h
