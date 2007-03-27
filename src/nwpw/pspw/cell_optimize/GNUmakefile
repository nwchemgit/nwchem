
# $Id: GNUmakefile,v 1.8 2007-03-27 02:02:48 bylaska Exp $
#
 OBJ_OPTIMIZE = \
		cell_optimize_input.o \
                cell_optimize.o \
                update_unitcell.o \
		pspw_stress_numerical.o \
		stress_fd.o \
		pspw_force_numerical.o \
		force_fd.o

	  OBJ =  

      LIBRARY = libpspw.a

    USES_BLAS = update_unitcell.F stress_fd.F force_fd
                
  LIB_DEFINES =


include ../../../config/makefile.h
include ../../../config/makelib.h

