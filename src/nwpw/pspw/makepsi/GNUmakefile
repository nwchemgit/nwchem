#
# $Id: GNUmakefile,v 1.11 2004-05-14 17:37:51 edo Exp $
#
 OBJ_OPTIMIZE = \
                wvfnc_init_input.o \
                wvfnc_init.o \
                wvfnc_new.o \
                wvfnc_expander_input.o \
                wvfnc_expander.o \
                wvfnc_reformat.o \
                v_wvfnc_init_input.o \
                v_wvfnc_init.o \
                v_wvfnc_new.o \
		pspw_wannier_input.o \
		pspw_wannier.o \
		silvestrelli_overlap.o \
                silvestrelli_minimize.o \
		expand_cell.o \
		wvfnc_expand_cell.o


	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = wvfnc_init.F v_wvfnc_init.F wvfnc_expander.F v_wvfnc_new.F wvfnc_new.F silvestrelli_minimize.F expand_cell.F  wvfnc_expand_cell.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

