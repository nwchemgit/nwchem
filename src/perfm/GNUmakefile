#
# $Id: GNUmakefile,v 1.3 2006-05-19 18:59:26 edo Exp $
#

  LIBRARY = libperfm.a
  
  OBJ = perfm_start.o perfm_end.o perfm_flop.o
  OBJ_OPTIMIZE = dist.o
  HEADERS = dist.fh

USES_BLAS = 

include ../config/makefile.h
include ../config/makelib.h
LIB_INCLUDES += -I$(MPI_INCLUDE)
 

testdist:	testdist.o $(LIBRARY_PATH)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ testdist.o $(LIBS)

