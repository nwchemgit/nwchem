#$Id: GNUmakefile,v 1.22 2003-09-15 01:24:25 sohirata Exp $
#
# CCSDTQ and CCSDTQ_LAMBDA will not be compiled by default.
# To include these, setenv CCSDTQ yes and/or setenv CCSDTQ_LAMBDA yes.

SUBDIRS = include ccd ccsd lccd lccsd qcisd ccsdt ccsdtq \
          cisd cisdt cisdtq mbpt1 mbpt2 mbpt3 mbpt4 \
          ccsd_lambda ccsdt_lambda ccsdtq_lambda \
          eomccsd

OBJ_OPTIMIZE = tce_energy.o tce_input.o tce_init.o tce_tile.o tce_ao1e.o tce_mo1e.o tce_ao2e.o \
      tce_filename.o createfile.o deletefile.o reconcilefile.o dratoga.o gatodra.o gatoeaf.o \
      get_block.o put_block.o add_block.o tce_tidy.o sf_test.o \
      tce_mo2e.o tce_mo1e_offset.o tce_mo2e_offset.o \
      tce_e_offset.o tce_t1_offset.o tce_t2_offset.o tce_t3_offset.o tce_t4_offset.o \
      tce_guess_t1.o tce_guess_t2.o tce_zero.o \
      tce_residual_t1.o tce_residual_t2.o tce_residual_t3.o tce_residual_t4.o \
      tce_jacobi_t1.o tce_jacobi_t2.o tce_jacobi_t3.o tce_jacobi_t4.o \
      tce_jacobi_l1.o tce_jacobi_l2.o tce_jacobi_l3.o tce_jacobi_l4.o \
      ma_pack.o sf_print.o ma_print_compact.o ma_transpose.o ma_scale.o \
      tce_diis_init.o tce_diis_tidy.o tce_diis_switch.o tce_diis.o \
      copyfile.o ddotfile.o daxfile.o daxpyfile.o dscalfile.o \
      tce_aod1.o tce_nud1.o tce_mod1.o tce_mod1_offset.o sym_abelian_axis.o \
      tce_eom_init.o tce_eom_tidy.o tce_eom_xguess.o tce_eom_xdiagon.o tce_eom_yguess.o tce_eom_ydiagon.o \
      tce_x1_offset.o tce_x2_offset.o tce_x3_offset.o tce_x4_offset.o \
      tce_y1_offset.o tce_y2_offset.o tce_y3_offset.o tce_y4_offset.o \
      tce_jacobi_x1.o tce_jacobi_x2.o tce_jacobi_x3.o tce_jacobi_x4.o \
      tce_jacobi_y1.o tce_jacobi_y2.o tce_jacobi_y3.o tce_jacobi_y4.o \
      parT_t3_fuse.o tce_mbpt2.o

LIB_INCLUDES = -I./include

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

