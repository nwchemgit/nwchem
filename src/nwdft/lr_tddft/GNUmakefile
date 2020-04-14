#$Id$

OBJ = tddft_energy.o tddft_input.o tddft_init.o tddft_memory.o tddft_iguess.o \
tddft_dump.o tddft_davidson.o tddft_transfm.o tddft_contract.o tddft_epsilon.o \
tddft_subspace.o tddft_diagon.o tddft_residual.o tddft_analysis.o tddft_dftinit.o \
tddft_finalize.o tddft_schmidt.o tddft_fold.o \
tddft_bdtensor.o tddft_patch.o tddft_civecs_phase.o tddft_oscstr.o \
tddft_checkpoint.o tddft_restart.o tddft_s2.o tddft_simpleso.o \
tddft_lanczos.o tddft_lanczos_tda.o tddft_dipole.o tddft_prod_v.o

OBJ_OPTIMIZE = 

LIBRARY = libnwdft.a

USES_BLAS = tddft_diagon.F tddft_transfm.F tddft_analysis.F tddft_restart.F tddft_civecs_phase.F

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

include ../../config/makefile.h
include ../../config/makelib.h

