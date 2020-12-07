OBJ = gw_input.o gw_defaults.o task_gw.o dft_energy_gw.o\
	gw_init.o gw_run.o gw_finalize.o gw_ri_init.o\
	gw_memory.o ga_chol_inv_seq.o gw_vxc.o\
	gw_cdgw.o

OBJ_OPTIMIZE =

LIBRARY = libnwdft.a

HEADERS = ./gw.fh ../grid/cgridfile.fh

LIB_DEFINES =
LIB_INCLUDES = -I../include -I../grid

USES_BLAS =

include ../../config/makefile.h
include ../../config/makelib.h
