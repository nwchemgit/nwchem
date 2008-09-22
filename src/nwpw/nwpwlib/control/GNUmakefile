#
# $Id: GNUmakefile,v 1.10 2008-09-22 22:55:15 bylaska Exp $
#

 OBJ_OPTIMIZE = \
            control.o \
			psi_get_ne.o \
			psi_get_version.o \
			psi_get_header.o \
			cpsi_get_header.o \
			get_psp_type.o


	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = control.F

  LIB_DEFINES = 

ifneq ($(USE_MPI),y)
  LIB_DEFINES += -DNOMPI
endif


include ../../../config/makefile.h
include ../../../config/makelib.h

