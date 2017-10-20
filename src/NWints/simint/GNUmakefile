# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../config/makefile.h

ifdef USE_SIMINT
      LIB_DEFINES += $(SIMINT_INCLUDE) -DUSE_SIMINT
      OBJ_OPTIMIZE = nwcsimint_test.o
else
      OBJ_OPTIMIZE = 
endif


  OBJ = nwcsim_facef90.o nwcsim_hf2.o
        
       HEADERS = 

       LIBRARY = libnwints.a

     USES_BLAS = 

include ../../config/makelib.h

#  ifeq ($(FC),efc)
#            FOPTIONS += -save
#           FOPTIMIZE += -O0 -g
#  endif

