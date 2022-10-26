#$Id$

 LIBRARY = libnwcutil.a


  OBJ_OPTIMIZE = geom_hnd.o

     OBJ = geom.o geom_input.o geom_input2.o geom_3d.o geom_2d.o geom_1d.o geom_numcore.o \
           geom_checksum.o geom_print_ecce.o  geom_freeze.o geom_fragment.o geom_getsym.o

 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F geom_print_ecce.F geom_hnd.F geom_3d.F geom_getsym.F geom_2d.F
 UNSET_OPENMP = 1
include ../config/makefile.h
include ../config/makelib.h


  ifeq ($(_FC),pgf90)
            FOPTIONS += -Mnorecursive
  endif
  ifeq ($(FC),gfortran)
        ifdef USE_OPENMP
            FOPTIONS += -fno-recursive -fno-openmp
	endif
  endif
  ifeq ($(_FC),gfortran)
	ifdef USE_OPENMP
            FOPTIONS += -fno-recursive -fno-openmp
	endif
  endif
  ifeq ($(FC),pathf90)
            FDEBUG += -O1 -g -LANG:heap_allocation_threshold=0
  endif
  ifeq ($(_CPU),ia64)
    ifeq ($(FC),ifort)
            FOPTIONS += -save
            FOPTIMIZE += -O0
    endif
  endif
  ifeq ($(_FC),ifort)
            FOPTIONS += -save
  endif
  ifeq ($(_FC),xlf)
	    FOPTIONS += -qsave
  endif
   ifeq ($(TARGET),HPUX64)
             FOPTIMIZE += +O0
   endif
   ifeq ($(TARGET),HPUX)
             FOPTIMIZE += +O0
   endif
