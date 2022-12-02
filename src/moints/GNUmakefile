#$Id$

include ../config/makefile.h
        OBJ = atwtest.o \
	      direct_mp2.o \
	      mointdriver.o \
	      moints_1e.o \
	      moints_dbg.o \
	      moints_gsh1.o \
	      moints_mem.o \
	      moints_mp2.o \
	      moints_stats.o \
              moints_full.o 

OBJ_OPTIMIZE = \
	      moints_cutil.o \
              moints_aodisk.o \
              moints_2x.o \
              moints_6x.o \
	      moints_semi.o \
              moints_aux2.o  

    LIBRARY = libmoints.a

  USES_BLAS = moints_2x.F \
              moints_6x.F \
              moints_aux2.F \
              moints_ccsd.F \
              moints_full.F \
              moints_semi.F 

LIB_DEFINES = -DDEBUG_PRINT

    HEADERS = moints_stats.fh
ifeq ($(ARMCI_NETWORK),MPI-PR)
   LIB_DEFINES += -DACC_STRIPS
endif

#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

include ../config/makelib.h

