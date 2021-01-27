OBJ = gw_input.o gw_defaults.o task_gw.o dft_energy_gw.o\
	gw_init.o gw_finalize.o \
	gw_memory.o ga_chol_inv_seq.o gw_vxc.o \
	gw_cdgw_newton.o gw_cdgw_graph.o gw_cdgw_linear.o

USES_BLAS = gw_cdgw.F gw_ri_init.F gw_run.F ga_chol_inv_seq.F\
            gw_lu_solve.F gw_invdiag ga_chol_inv_seq.F\
	    gw_buildrn.o gw_buildiw.o

OBJ_OPTIMIZE = gw_ri_init.o gw_cdgw.o gw_run.o gw_lu_solve.o\
               gw_invdiag.o gw_buildrn.o gw_buildin.o\
	       gw_buildiw.o

LIBRARY = libnwdft.a

HEADERS = 

LIB_DEFINES =
LIB_INCLUDES = -I../include -I../grid -I./

include ../../config/makefile.h
include ../../config/makelib.h
