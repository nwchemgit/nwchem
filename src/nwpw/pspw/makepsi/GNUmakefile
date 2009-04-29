#
# $Id: GNUmakefile,v 1.13 2007-01-02 18:36:52 bylaska Exp $
#
 OBJ_OPTIMIZE = \
                wvfnc_init_input.o \
                wvfnc_init.o \
                wvfnc_new.o \
                wvfnc_expander_input.o \
                wvfnc_expander.o \
                wvfnc_reformat.o \
                wvfnc_ionize.o \
                wvfnc_adjust.o \
                wvfnc_convert_restricted.o \
                v_wvfnc_init_input.o \
                v_wvfnc_init.o \
                v_wvfnc_new.o \
		pspw_wannier_input.o \
		pspw_wannier.o \
		silvestrelli_overlap.o \
                silvestrelli_minimize.o \
		expand_cell.o \
		wvfnc_expand_cell.o \
                band_reformat_c_wvfnc.o


	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = wvfnc_init.F v_wvfnc_init.F wvfnc_expander.F v_wvfnc_new.F wvfnc_new.F silvestrelli_minimize.F expand_cell.F  wvfnc_expand_cell.F band_reformat_c_wvfnc.F wvfnc_adjust.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

