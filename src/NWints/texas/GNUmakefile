# $Id: GNUmakefile,v 1.8 1995-11-01 17:16:18 d3e129 Exp $
 
ifeq ($(USE_TEXAS),YEP)
      LIB_DEFINES = -DUSE_TEXAS
  OBJ_OPTIMIZE = texas_face.o\
                 texas_hf.o\
                 big_service.o\
                 getmem.o\
                 amshift.o\
                 assemblx.o\
            spec_blocking.o\
            spec_calcint.o\
                 cshneg.o\
               d_datlog.o\
            spec_destiny.o\
                 memoha.o\
                 pak.o\
                 precalc2.o\
                 prepint.o\
           label_reorder.o\
                 rt123.o\
                 shells.o\
                 time.o\
                 trobsa.o\
                 zeroint.o \
		 missing.o
else
      OBJ = texas_stubs.o
      OBJ_OPTIMIZE = 
endif

       LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h
