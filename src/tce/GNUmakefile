#$Id: GNUmakefile,v 1.14 2003-01-31 04:36:42 sohirata Exp $

OBJ_OPTIMIZE = tce_energy.o tce_input.o tce_init.o tce_tile.o tce_ao1e.o tce_mo1e.o tce_ao2e.o \
      tce_filename.o createfile.o deletefile.o reconcilefile.o \
      get_block.o put_block.o add_block.o tce_tidy.o sf_test.o \
      tce_mo2e.o tce_mo1e_offset.o tce_mo2e_offset.o \
      tce_e_offset.o tce_t1_offset.o tce_t2_offset.o tce_t3_offset.o tce_t4_offset.o \
      tce_guess_t1.o tce_guess_t2.o tce_guess_t3.o tce_guess_t4.o tce_zero.o \
      tce_residual_t1.o tce_residual_t2.o tce_residual_t3.o tce_residual_t4.o \
      tce_jacobi_t1.o tce_jacobi_t2.o tce_jacobi_t3.o tce_jacobi_t4.o \
      ma_pack.o sf_print.o ma_print_compact.o ma_transpose.o ma_scale.o \
      ccd_e.o ccd_t2.o lccd_e.o lccd_t2.o ccsd_e.o ccsd_t1.o ccsd_t2.o lccsd_e.o lccsd_t1.o lccsd_t2.o \
      ccsdt_e.o ccsdt_t1.o ccsdt_t2.o ccsdt_t3.o ccsdtq_e.o ccsdtq_t1.o ccsdtq_t2.o ccsdtq_t3.o ccsdtq_t4.o \
      cisd_e.o cisd_c1.o cisd_c2.o cisdt_e.o cisdt_c1.o cisdt_c2.o cisdt_c3.o \
      cisdtq_e.o cisdtq_c1.o cisdtq_c2.o cisdtq_c3.o cisdtq_c4.o \
      tce_mbpt2.o mbpt1_t1.o mbpt1_t2.o mbpt2_e.o mbpt2_t1.o mbpt2_t2.o mbpt2_t2_b.o mbpt2_t3.o mbpt2_t4.o \
      mbpt3_e.o mbpt3_t1.o mbpt3_t2.o mbpt4_e.o qcisd_e.o qcisd_t1.o qcisd_t2.o \
      tce_diis_init.o tce_diis_tidy.o tce_diis_switch.o tce_diis_t1.o tce_diis_t2.o tce_diis_t3.o tce_diis_t4.o \
      copyfile.o ddotfile.o daxpyfile.o

OBJ = 

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

