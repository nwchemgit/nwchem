#$Id: GNUmakefile,v 1.32 1995-06-16 20:51:14 rg240 Exp $

        OBJ = atwtest.o \
	      direct_mp2.o \
	      ga_create_JKblocked.o \
	      ga_check_JKblocked.o \
	      mointdriver.o \
	      moints_dbg.o \
	      moints_gsh1.o \
	      moints_mem.o \
	      moints_mp2.o \
	      moints_1e.o \
	      moints_pascii.o \
	      moints_stats.o 

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
              moints_aux2.F direct_mp2.F mointdriver.F

LIB_DEFINES =  -DDEBUG_PRINT

#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#


include ../config/makefile.h
include ../config/makelib.h

#
#
#
moints_2x.o:	moints_stats.fh
moints_6x.o:	moints_stats.fh
moints_stats.o:	moints_stats.fh

















