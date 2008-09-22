#
# $Id: GNUmakefile,v 1.9 2008-09-22 22:55:16 bylaska Exp $
#
 SUBDIRS = pspw_LJ pspw_Q 
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

