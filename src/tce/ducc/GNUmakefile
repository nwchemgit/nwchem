#$Id$

OBJ_OPTIMIZE =  map_hcore.o rot_vir.o ducc.o hh.o hp.o pp.o hhhh.o pppp.o hhhp.o ppph.o hhpp.o hphp.o
include ../../config/makefile.h

#ifeq ($(_FC),ifort)
#FOPTIMIZE = -O3 -ip -unroll
#endif

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = rot_vir.F

LIB_DEFINES = -DDEBUG_PRINT

ifdef OPTIMIZE_DUCC
LIB_DEFINES += -DOPTIMIZE_DUCC
endif
#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

HEADERS = 

include ../../config/makelib.h

