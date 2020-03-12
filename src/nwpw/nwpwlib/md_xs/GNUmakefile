#
# $Id: GNUmakefile 25011 2013-12-19 17:45:25Z d3y133 $
#

ifndef USE_FEFF
  SUBDIRS = FEFF6L
endif

 OBJ_OPTIMIZE = \
		md_xs_input.o \
		md_xs.o

	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = md_xs.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h


ifeq ($(USE_FEFF),y)
  LIB_DEFINES += -DFEFF
  ifeq ($(USE_FEFF8),y)
     LIB_DEFINES += -DFEFF8
  endif
else
  LIB_DEFINES += -DFEFF
  LIB_DEFINES += -DFEFF8
endif

