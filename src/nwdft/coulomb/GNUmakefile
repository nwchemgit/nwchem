# $Id: GNUmakefile,v 1.27 1997-07-21 20:23:15 d3h449 Exp $
ifeq ($(TARGET),SGITFP)
      OBJ_OPTIMIZE = dft_getvc.o dft_n3c.o
      OBJ =  dft_3cincor.o\
	dft_fitcd.o				   dft_fitvc.o	\
	dft_get2eri.o		dft_genspm.o dft_getspm.o 
else
#      OBJ = dft_genspm.o dft_getspm.o 
#      OBJ_OPTIMIZE = dft_fitcd.o dft_3cincor.o dft_n3c.o \
#	dft_fitvc.o	\
#	dft_get2eri.o			   dft_getvc.o	
#

	OBJ_OPTIMIZE = \
	dft_3cincor.o \
	dft_fitcd.o \
	dft_fitvc.o \
	dft_genspm.o \
	dft_get2eri.o \
	dft_getspm.o \
	dft_getvc.o \
	dft_n3c.o


endif

      LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

 LIB_INCLUDES = -I../include



      USES_BLAS = dft_fitcd.F dft_fitvc.F dft_genspm.F dft_3cincor.F

include ../../config/makefile.h
include ../../config/makelib.h

