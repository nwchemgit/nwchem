
        OBJ = atwtest.o \
              brdcst_vec.o \
	      ga_create_moblocked.o \
	      ga_check_moblocked.o \
	      mointdriver.o \
	      moints_mp2.o \
	      moints_xx.o \
	      moints.o 

    LIBRARY = libmoints.a

  USES_BLAS = moints.F moints_xx.F

 LIB_DEFINES =  -DDEBUG_PRINT -DDEBUG_PARALLEL

#	       -DINTS_ONLY
#              -DMOINTS_INTS_DUMMY
#              -DDEBUG_PRINT
#              -DDEBUG_PARALLEL
#	       -DPERMUT_SYMM
#              -DOLD_ALGORITHM


include ../config/makefile.h
include ../config/makelib.h

















