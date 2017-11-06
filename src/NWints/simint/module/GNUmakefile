# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../../config/makefile.h

  OBJ_OPTIMIZE = nwcsim_module.o
ifdef USE_SIMINT
      LIB_DEFINES += -I$(SIMINT_HOME)/include/simint -DUSE_SIMINT
endif

ifdef SIMINT_GRADIENT
      LIB_DEFINES +=  -DSIMINT_GRADIENT
endif

        
       HEADERS = 

       LIBRARY = libnwints.a

     USES_BLAS = 

include ../../../config/makelib.h


