#$Id$

OBJ          = cr_ccsd_t.o lr_ccsd_t.o 

OBJ_OPTIMIZE = ccsd_t.o ccsd_t_doubles.o ccsd_t_singles.o lambda_ccsd_t.o lambda_ccsd_t_left.o\
               ccsd_t_doubles_l.o ccsd_t_singles_l.o \
               ccsd2_t.o ccsd2_t_left.o ccsd2_q.o ccsd2_q_left.o ccsd2_q_right.o \
               ccsdt2_q.o ccsdt2_q_left.o ccsdt2_q_right.o \
               cr_ccsd_t_D.o cr_ccsd_t_E.o cr_ccsd_t_N.o lr_ccsd_tq.o \
               cr_qcisd_t.o m2_dot_m1t1.o ccsd_t_restart.o  \
               ccsd_t_6dts.o ccsd_t_singles_l_6dts.o ccsd_t_doubles_l_6dts.o \
	           ccsd_t_dts_kernels_omp.o ccsd_t_transpose.o \
               offl_ccsd_t_doubles_l.o ccsd_t_kernels_omp.o \
               tce_hashnsort.o ccsd_t_pstat.o \
               ccsd_t_dot.o ccsd_t_neword.o
ifdef USE_OFFLOAD
  OBJ_OPTIMIZE += offl_ccsd_t_singles_l.o offl_ccsd_t_singles_l_1.o
endif

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_t_doubles.F ccsd_t_singles.F \
            ccsd_t_doubles_l.F ccsd_t_singles_l.F \
            lambda_ccsd_t_left.F \
            ccsd2_q_left.F ccsd2_q_right.F ccsd2_t_left.F \
            ccsdt2_q_left.F ccsdt2_q_right.F \
            cr_ccsd_t_D.F cr_ccsd_t_E.F cr_ccsd_t_N.F lr_ccsd_tq.F \
            cr_qcisd_t.F m2_dot_m1t1.F \
            ccsd_t_singles_l_6dts.F ccsd_t_doubles_l_6dts.F ccsd_t.F


ifdef TCE_CUDA
 LIB_INCLUDES += -I../include 
 OBJ_OPTIMIZE += hybrid.o memory.o ccsd_t_gpu.o ccsd_t_singles_gpu.o ccsd_t_doubles_gpu.o
 ifdef USE_TTLG
  LIB_INCLUDES += -I../ttlg/includes
  OBJ_OPTIMIZE +=  sd_t_total_ttlg.o
 else
  OBJ_OPTIMIZE +=  sd_t_total.o
 endif
 USES_BLAS += ccsd_t_singles_gpu.F ccsd_t_doubles_gpu.F
 LIB_DEFINES += $(CUDA_INCLUDE)
endif

ifdef TCE_HIP
 LIB_INCLUDES = -I../include
 OBJ_OPTIMIZE += hybrid.o memory.o ccsd_t_gpu.o ccsd_t_singles_gpu.o ccsd_t_doubles_gpu.o sd_t_total.o
 USES_BLAS += ccsd_t_singles_gpu.F ccsd_t_doubles_gpu.F
 LIB_DEFINES += $(HIP_INCLUDE)
endif


ifdef TCE_OPENACC
  OBJ_OPTIMIZE += ccsd_t_kernels_acc.o
endif

ifdef TCE_NEW_OPENMP
 OBJ_OPTIMIZE += ccsd_t_omp.o
endif

LIB_DEFINES += -DDEBUG_PRINT

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

ifdef USE_OFFLOAD
FOPTIMIZE += $(OFFLOAD_FOPTIONS)
endif
