
# $Id$
#
include ../../../config/makefile.h
 OBJ_OPTIMIZE = \
		cell_optimize_input.o \
                cell_optimize.o \
                update_unitcell.o \
		pspw_stress_numerical.o \
		pspw_force_numerical.o \
		force_fd.o

  ifeq ($(_IFCV10),Y)
          OBJ = stress_fd.o
  else
          OBJ_OPTIMIZE += stress_fd.o
  endif


      LIBRARY = libpspw.a

    USES_BLAS = update_unitcell.F stress_fd.F force_fd.F
                
  LIB_DEFINES =


include ../../../config/makelib.h

