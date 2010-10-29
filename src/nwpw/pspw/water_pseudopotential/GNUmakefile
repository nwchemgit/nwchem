
# $Id$
#
# intf and aux break on IBM SP xlhpf with optimization ... why? 
 OBJ_OPTIMIZE = \
                water_pseudopotential_input.o \
				waterpsp.o \
				LJparam_Ion.o \
				LJ_Ion.o \
				BLCJ_SR.o \
				BLCJ_LR.o \
				BLCJ_Coulomb.o \
				BLCJ_Polarization.o \
				BLCJ_Intra.o \
				BLCJ_Inter.o

	  OBJ = 


      LIBRARY = libpspw.a

    USES_BLAS = waterpsp.F BLCJ_SR.F
                
  LIB_DEFINES =


include ../../../config/makefile.h
include ../../../config/makelib.h

