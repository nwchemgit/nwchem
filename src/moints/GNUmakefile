#$Id: GNUmakefile,v 1.28 1995-02-02 23:18:43 d3g681 Exp $

        OBJ = atwtest.o \
              brdcst_vec.o \
	      direct_mp2.o \
	      ga_create_JKblocked.o \
	      ga_check_JKblocked.o \
	      mointdriver.o \
	      moints_dbg.o \
	      moints_gsh1.o \
	      moints_mem.o \
	      moints_mp2.o \
	      moints_stats.o 

OBJ_OPTIMIZE = \
	      moints_2x.o \
              moints_6x.o \
              moints_aux2.o 

    LIBRARY = libmoints.a

  USES_BLAS = \
	      moints_2x.F \
              moints_6x.F \
              moints_aux2.F 

LIB_DEFINES =  -DDEBUG_PRINT

#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#              -DMOINTS_STATS         /* Include statistics gathering */
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

















