#$Id: GNUmakefile,v 1.50 2004-01-02 21:38:26 edo Exp $

 LIBRARY = libnwcutil.a


  OBJ_OPTIMIZE = geom_hnd.o

     OBJ = geom.o geom_input.o geom_input2.o geom_3d.o geom_2d.o geom_1d.o geom_numcore.o \
           geom_checksum.o geom_print_ecce.o  geom_freeze.o geom_symmol.o

 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F geom_print_ecce.F geom_hnd.F
include ../config/makefile.h
include ../config/makelib.h


  ifeq ($(FC),efc)
	    FOPTIONS += -save
            FOPTIMIZE += -O0
  endif
  ifeq ($(FC),ifc)
	    FOPTIONS += -save
  endif
  ifeq ($(FC),xlf)
	    FOPTIONS += -qsave
  endif
