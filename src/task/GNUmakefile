#
# $Id$
#

  LIBRARY = libnwctask.a
  
  OBJ = task_energy.o task_gradient.o task_optimize.o task_freq.o \
	task_hessian.o task_dynamics.o task_thermo.o task_et.o\
        task_save_state.o task.o task_input.o task_shell_input.o task_num_grad.o \
        task_vscf.o task_property.o task_dntmc.o task_bsse.o \
        task_jefftce.o task_ncc.o task_rfcf.o \
	task_gwmol.o \
	task_aoints.o task_nadtrc.o

USES_BLAS = task_num_grad.F task_bsse.F

ifdef SKIP_ECCE
 LIB_DEFINES += -DSKIP_ECCE
endif

include ../config/makefile.h
include ../config/makelib.h
