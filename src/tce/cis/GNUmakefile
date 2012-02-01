#$Id$

OBJ_OPTIMIZE =  cis.o cis_x1.o cis_xdiagon.o cis_xtidy.o

# For unknown reasons the IBM xlf compiler generated an ICE on Power7 systems
# while trying to optimize this code. This problem might go away in future.
# Date 01/11/2012.
ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),LAPI64 IBM64))
  OBJ += cis_xguess.o
else
  OBJ_OPTIMIZE += cis_xguess.o
endif

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cis_x1.F cis_xdiagon.F \


LIB_DEFINES = -DDEBUG_PRINT
#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

HEADERS = 

include ../../config/makefile.h
include ../../config/makelib.h

