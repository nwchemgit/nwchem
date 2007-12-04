
# $Id: GNUmakefile,v 1.9 2007-12-04 18:21:52 d3p307 Exp $
#
 OBJ_OPTIMIZE = \
		cell_optimize_input.o \
                cell_optimize.o \
                update_unitcell.o \
		pspw_stress_numerical.o \
		pspw_force_numerical.o \
		force_fd.o

ifeq ($(NWCHEM_TARGET),LINUX)
  ifeq ($(FC),ifort)
          OBJ = stress_fd.o
  else
          OBJ_OPTIMIZE += stress_fd.o
  endif
else
          OBJ_OPTIMIZE += stress_fd.o
endif

      LIBRARY = libpspw.a

    USES_BLAS = update_unitcell.F stress_fd.F force_fd
                
  LIB_DEFINES =


include ../../../config/makefile.h
include ../../../config/makelib.h

