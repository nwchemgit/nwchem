# $Id: GNUmakefile,v 1.32 2004-01-20 19:08:42 edo Exp $
ifeq ($(TARGET),SGITFP)
# BEGIN_WIN32_IGNORE
      OBJ_OPTIMIZE = dft_getvc.o dft_n3c.o
      OBJ =  dft_3cincor.o dft_mem3c.o\
	dft_fitcd.o				   dft_fitvc.o	\
	dft_get2eri.o		dft_genspm.o dft_getspm.o 
# END_WIN32_IGNORE
else
#      OBJ = dft_genspm.o dft_getspm.o 
#      OBJ_OPTIMIZE = dft_fitcd.o dft_3cincor.o dft_n3c.o \
#	dft_fitvc.o	\
#	dft_get2eri.o			   dft_getvc.o	
#

	OBJ= \
	dft_3cincor.o  dft_mem3c.o\
	dft_fitcd.o \
	dft_fitvc.o \
	dft_genspm.o \
	dft_get2eri.o \
	dft_getspm.o \
	dft_getvc.o \
	dft_n3c.o


endif

      LIBRARY = libnwdft.a

 LIB_DEFINES =

 LIB_INCLUDES = -I../include



      USES_BLAS = dft_fitcd.F dft_fitvc.F dft_genspm.F dft_3cincor.F dft_get2eri.F

include ../../config/makefile.h
include ../../config/makelib.h

