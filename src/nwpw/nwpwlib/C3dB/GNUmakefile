#
# $Id: GNUmakefile,v 1.3 2003-10-09 00:36:06 edo Exp $
#

 OBJ_OPTIMIZE = C3dB-new.o
	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = C3dB-new.F

  LIB_DEFINES = 


include ../../../config/makefile.h
ifeq ($(FC),pgf90)
 OBJ = C3dB-new.o
 OBJ_OPTIMIZE = 
endif
include ../../../config/makelib.h

