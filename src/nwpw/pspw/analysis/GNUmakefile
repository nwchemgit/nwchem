#
# $Id: GNUmakefile,v 1.6 2004-05-14 16:43:58 edo Exp $
#
     SUBDIRS = APC


 OBJ_OPTIMIZE = \
		analysis_input.o \
		analysis_expansion_coef.o \
		analysis.o \
                orb_pop.o \
                orb_pop_kawai.o \
                orb_pop_aorb.o \
                orb.o \
                pspw_Lubin_water_analysis.o 


      LIBRARY = libpspw.a

    USES_BLAS =  pspw_Lubin_water_analysis.F analysis.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

