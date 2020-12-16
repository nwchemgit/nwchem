OBJ = gw_input.o gw_defaults.o task_gw.o dft_energy_gw.o\
	gw_init.o gw_finalize.o \
	gw_memory.o ga_chol_inv_seq.o gw_vxc.o

USES_BLAS = gw_cdgw.F gw_ri_init.F gw_run.F ga_chol_inv_seq.F

OBJ_OPTIMIZE = gw_ri_init.o gw_cdgw.o gw_run.o

LIBRARY = libnwdft.a

HEADERS = ./gw.fh ../grid/cgridfile.fh

LIB_DEFINES =
LIB_INCLUDES = -I../include -I../grid

USES_BLAS =

include ../../config/makefile.h
include ../../config/makelib.h
