#$Id: GNUmakefile,v 1.45 1996-01-26 21:59:47 d3g681 Exp $

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
              moints_full.o

ifneq ($(TARGET),CRAY-T3D)
        OBJ += crayio.o
endif

OBJ_OPTIMIZE = \
	      moints_2x.o \
              moints_6x.o \
	      moints_ccsd.o \
              moints_aux2.o 

    LIBRARY = libmoints.a

  USES_BLAS = \
	      moints_2x.F \
              moints_6x.F \
	      moints_ccsd.F \
              moints_aux2.F \
              direct_mp2.F \
              moints_semi.F \
              mointdriver.F

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

