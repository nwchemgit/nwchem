#
# $Id: GNUmakefile,v 1.10 2009-03-26 00:19:01 bylaska Exp $
#
 SUBDIRS = pspw_LJ pspw_Q mmq
 OBJ_OPTIMIZE = \
                pspw_qmmm_input.o \
                pspw_qmmm.o \
                pspw_2qmmm.o \
		pspw_springs.o

      OBJ =


      LIBRARY = libpspw.a

    USES_BLAS = pspw_qmmm.F pspw_qmmm2.F
  LIB_DEFINES =
  LIB_INCLUDES = 

include ../../../config/makefile.h
include ../../../config/makelib.h

