#
# $Id: GNUmakefile,v 1.4 2001-11-28 22:15:22 bylaska Exp $
#
# intf and aux break on IBM SP xlhpf with optimization ... why? 
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
				getfilling.o


	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = wvfnc_init.F v_wvfnc_init.F wvfnc_expander.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

