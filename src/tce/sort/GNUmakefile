#$Id: GNUmakefile 26511 2014-12-10 17:54:10Z jhammond $
#

OBJ_OPTIMIZE = tce_sort_4kg.o

# OLD DEFAULT
#OBJ_OPTIMIZE += tce_sort.o tce_sortacc.o

# OLD NON-DEFAULT (missing some cases that can be added per-file as below)
#OBJ_OPTIMIZE += tce_sortacc_saday.o tce_sort_saday.o

# These are pointless utility routines that exist only because of code generator.
OBJ_OPTIMIZE += tce_sort0.o tce_sortacc0.o

# NEW sorts - faster, simpler, amenable to OpenMP
OBJ_OPTIMIZE += new_sort2.o 
OBJ_OPTIMIZE += new_sort4.o 
OBJ_OPTIMIZE += new_sort6.o 
# OLD sorts
#OBJ_OPTIMIZE += tce_sort2.o tce_sortacc2.o
#OBJ_OPTIMIZE += tce_sort4.o tce_sortacc4.o
#OBJ_OPTIMIZE += tce_sort6.o tce_sortacc6.o

# These have not been optimized yet because quadruples do not get used much.
OBJ_OPTIMIZE += tce_sort8.o tce_sortacc8.o

# This is used in the EOM code.
OBJ_OPTIMIZE += tce_sort1.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = 

HEADERS =

# This has to go before makefile.h
#OPTIMIZE = 1

include ../../config/makefile.h

# This has to go after makefile.h
#FOPTIMIZE += -qopt-streaming-stores=always

include ../../config/makelib.h

