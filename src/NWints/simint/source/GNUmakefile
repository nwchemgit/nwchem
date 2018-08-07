# $Id: GNUmakefile 26876 2015-02-24 06:32:05Z edo $

include ../../../config/makefile.h

  OBJ_OPTIMIZE = nwcsim_facef90.o nwcsim_hf2.o
ifdef USE_SIMINT
      LIB_DEFINES += -I$(SIMINT_HOME)/include/simint -DUSE_SIMINT
      LIB_DEFINES += -I../module
      OBJ_OPTIMIZE += nwcsimint_test.o

      SIMINT_OSTEI_MAXAM=$(shell grep SIMINT_OSTEI_MAXAM $(SIMINT_HOME)/include/simint/ostei/ostei_config.h|head -1|cut -c 28)

      SIMINT_OSTEI_MAXDER=$(shell grep SIMINT_OSTEI_MAXDER $(SIMINT_HOME)/include/simint/ostei/ostei_config.h|head -1|cut -c 29)

      SIMINT_MAXDER_GT_0 = $(shell [ $(SIMINT_OSTEI_MAXDER) -gt 0  ] && echo true)

      LIB_DEFINES +=  -DSIMINT_MAXAM=$(SIMINT_OSTEI_MAXAM)
      ifeq ($(SIMINT_MAXDER_GT_0),true)
           LIB_DEFINES +=  -DSIMINT_GRADIENT
      endif
endif


        
       HEADERS = 

       LIBRARY = libnwints.a

     USES_BLAS = 

include ../../../config/makelib.h
ifdef USE_SIMINT
$(LIBRARY_PATH)(nwcsim_hf2.o):	../module/nwcsim_shell.mod
endif

#  ifeq ($(FC),efc)
#            FOPTIONS += -save
#           FOPTIMIZE += -O0 -g
#  endif

