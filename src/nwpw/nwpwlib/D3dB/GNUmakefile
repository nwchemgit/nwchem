#
# $Id: GNUmakefile,v 1.3 2003-10-09 00:36:06 edo Exp $
#

 OBJ_OPTIMIZE = D3dB-new.o
 OBJ = 


      LIBRARY = libnwpwlib.a

    USES_BLAS = D3dB-new.F

  LIB_DEFINES = 



include ../../../config/makefile.h
ifeq ($(FC),pgf90)
 OBJ_OPTIMIZE = 
 OBJ = D3dB-new.o
endif
include ../../../config/makelib.h

