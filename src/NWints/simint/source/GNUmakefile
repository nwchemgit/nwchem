# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../../config/makefile.h

  OBJ_OPTIMIZE = nwcsim_facef90.o nwcsim_hf2.o
ifdef USE_SIMINT
      LIB_DEFINES += -I$(SIMINT_HOME)/include/simint -DUSE_SIMINT
      LIB_DEFINES += -I../module
      OBJ_OPTIMIZE += nwcsimint_test.o
endif

ifdef SIMINT_GRADIENT
      LIB_DEFINES +=  -DSIMINT_GRADIENT
endif

        
       HEADERS = 

       LIBRARY = libnwints.a

     USES_BLAS = 

include ../../../config/makelib.h

#  ifeq ($(FC),efc)
#            FOPTIONS += -save
#           FOPTIMIZE += -O0 -g
#  endif

