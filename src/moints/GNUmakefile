#$Id: GNUmakefile,v 1.62 1997-04-01 02:04:58 rg240 Exp $

        OBJ = atwtest.o \
	      direct_mp2.o \
	      mointdriver.o \
	      moints_1e.o \
	      moints_dbg.o \
	      moints_gsh1.o \
	      moints_mem.o \
	      moints_mp2.o \
	      moints_semi.o \
	      moints_stats.o \
              moints_full.o \
              sifs.o

#              twopdm_backt.o 


ifneq ($(TARGET),CRAY-T3D)
        OBJ += crayio.o
endif

OBJ_OPTIMIZE = \
              moints_aodisk.o \
	      moints_cutil.o \
	      moints_2x.o \
              moints_6x.o \
              moints_aux2.o 

    LIBRARY = libmoints.a

USES_BLAS = moints_2x.F moints_6x.F moints_aux2.F moints_ccsd.F moints_full.F moints_semi.F sifs.F

LIB_DEFINES =  -DDEBUG_PRINT

    HEADERS = moints_stats.fh

#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#

include ../config/makefile.h
include ../config/makelib.h

