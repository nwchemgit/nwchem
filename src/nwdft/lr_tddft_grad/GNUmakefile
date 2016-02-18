#$Id$

OBJ = tddft_gradient.o tddft_grad_init.o tddft_grad_input.o tddft_grad_finalize.o \
      tddft_grad_print_info.o tddft_grad_ints.o tddft_grad_print_header.o \
      tddft_grad_e.o tddft_grad_root_e.o tddft_grad_mo.o tddft_grad_xy.o\
      tddft_grad_load_mo.o tddft_grad_load_xy.o tddft_grad_tp.o tddft_grad_compute_t.o \
      tddft_grad_r.o tddft_grad_compute_r.o \
      tddft_grad_trans_mo2ao.o tddft_grad_trans_ao2mo.o \
      tddft_grad_mat_to_vec.o tddft_grad_solve.o tddft_grad_epsilon.o \
      tddft_grad_solve_product.o tddft_grad_solve_precond.o \
      tddft_grad_compute_p.o tddft_grad_compute_w.o \
      tddft_grad_g.o \
      tddft_grad_compute_dao.o tddft_grad_get_dens.o tddft_grad_compute_g.o \
      tddft_grad_grad_1e.o tddft_grad_grad_2e.o \
      tddft_grad_make_2pdm.o tddft_grad_print_result.o tddft_grad_print_array.o \
      tddft_grad_compute_energy.o tddft_grad_anti_symmetrize.o \
      tddft_grad_gxc.o tddft_grad_irreg_acc.o tddft_grad_fxc_l.o \
      tddft_grad_fix_exc.o tddft_grad_irreg_acc.o \
      tddft_grad_store_gradient.o tddft_grad_util_print.o \
      tddft_grad_solve_load_x.o tddft_grad_solve_save_x.o \
      tddft_grad_dxcdmat.o tddft_grad_dfxc.o tddft_grad_quadv0b.o \
      tddft_grad_dvxc.o tddft_grad_mpole.o

OBJ_OPTIMIZE = 

LIBRARY = libnwdft.a

USES_BLAS = tddft_grad_gxc.F tddft_grad_compute_w.F tddft_grad_dxcdmat.F \
            tddft_grad_dfxc.F tddft_grad_quadv0b.F tddft_grad_dvxc.F tddft_grad_compute_g.F

LIB_DEFINES  = -DDEBUG_PRINT
LIB_INCLUDES = -I../include -I../grid 
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
