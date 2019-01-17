# $Id$
  LIBRARY = lib64to32.a
  
  OBJ_OPTIMIZE = yscal.o ycopy.o ydot.o ygemm.o ygemv.o \
                 iyamax.o ygesv.o ygesvd.o yaxpy.o yspsvx.o \
                 yposv.o ysyev.o yger.o yrot.o ylaset.o yswap.o\
                 ygebak.o ygebal.o ynrm2.o ylacpy.o ylascl.o\
                 ylange.o ytrevc.o yhseqr.o yorghr.o ygehrd.o\
                 ysterf.o ylarnv.o ylagtf.o ylagts.o yasum.o \
                 ypotri.o ypotrf.o ysygv.o ygeev.o ygeevx.o \
                 ifily.o\
                 jscal.o jaxpy.o jgemm.o jheev.o jcopy.o jdotc.o \
                 ijamax.o jgesvd.o

        ifndef BLAS_SIZE
	LIB_DEFINES += -DUSE_INTEGER8
	endif
        ifeq ($(BLAS_SIZE),8)
	LIB_DEFINES += -DUSE_INTEGER8
	endif

include ../config/makefile.h
include ../config/makelib.h

