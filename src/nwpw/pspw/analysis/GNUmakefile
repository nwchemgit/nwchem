#
# $Id$
#
     SUBDIRS = APC


 OBJ_OPTIMIZE = \
		analysis_input.o \
		analysis.o \
                orb_pop.o \
                orb_pop_kawai.o \
                orb_pop_aorb.o \
                orb.o \
                pspw_Lubin_water_analysis.o \
                pspw_atom_analysis.o 


      LIBRARY = libpspw.a

    USES_BLAS =  pspw_Lubin_water_analysis.F pspw_atom_analysis.F analysis.F orb_pop_kawai.F orb_pop_aorb.F


  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

