
        OBJ = atwtest.o \
              brdcst_vec.o \
	      ga_create_JKblocked.o \
	      ga_check_JKblocked.o \
	      mointdriver.o \
	      moints.o \
	      moints_dbg.o \
	      moints_shgrp.o \
	      moints_mp2.o


OBJ_OPTIMIZE =  moints_auxil.o 

    LIBRARY = libmoints.a

  USES_BLAS = moints.F 

LIB_DEFINES =  -DDEBUG_PRINT

#
# Possible #defines
#
#              -DDEBUG_PRINT
#              -DMOINTS_STATS
#	       -DMOINTS_PAR_STATS
#	       -DUSEGEMM
#


include ../config/makefile.h
include ../config/makelib.h

















