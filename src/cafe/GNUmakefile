#
OBJ_OPTIMIZE=cf_d.o cf_f.o\
             cf_finish.o\
             cf_forces.o\
             cf_init.o\
             cf_list.o\
             cf_pbc.o\
             cf_rdtop.o

LIBRARY=libnwargos.a

HEADER=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h