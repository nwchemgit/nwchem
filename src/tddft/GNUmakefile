#$Id: GNUmakefile,v 1.2 2002-04-12 22:55:51 sohirata Exp $

OBJ = tddft_energy.o tddft_input.o tddft_init.o tddft_memory.o tddft_iguess.o \
tddft_dump.o tddft_davidson.o tddft_transfm.o tddft_contract.o tddft_epsilon.o \
tddft_subspace.o tddft_diagon.o tddft_residual.o tddft_analysis.o tddft_dftinit.o \
tddft_finalize.o tddft_sort.o tddft_schmitt.o tddft_dgeev.o \
ma_print.o sym_irrepname.o sym_product.o

OBJ_OPTIMIZE = 

LIBRARY = libtddft.a

USES_BLAS = 

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

include ../config/makefile.h
include ../config/makelib.h

