#$Id: GNUmakefile,v 1.4 2002-10-22 16:23:02 sohirata Exp $

OBJ = tce_energy.o tce_input.o tce_init.o tce_tile.o tce_ao1e.o tce_mo1e.o tce_ao2e.o \
      tce_mo2e.o tce_mo1e_offset.o tce_mo2e_offset.o tce_e_offset.o tce_t1_offset.o tce_t2_offset.o tce_t3_offset.o \
      createfile.o deletefile.o tce_guess_t1.o tce_guess_t2.o tce_guess_t3.o tce_mbpt2.o tce_zero.o \
      tce_residual_t1.o tce_residual_t2.o tce_residual_t3.o tce_jacobi_t1.o tce_jacobi_t2.o tce_jacobi_t3.o \
      get_block.o put_block.o add_block.o tce_tidy.o \
      ma_pack.o sf_print.o ma_print_compact.o ma_transpose.o ma_scale.o \
      ccd_e.o ccd_t2.o lccd_e.o lccd_t2.o ccsd_e.o ccsd_t1.o ccsd_t2.o lccsd_e.o lccsd_t1.o lccsd_t2.o \
      ccsdt_e.o ccsdt_t1.o ccsdt_t2.o ccsdt_t3.o 

OBJ_OPTIMIZE = 

LIBRARY = libtce.a

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

