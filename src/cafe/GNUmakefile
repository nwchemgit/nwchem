#
# $Id$
#
OBJ_OPTIMIZE=cf_d.o cf_fp.o cf_f.o\
             cf_finish.o\
             cf_forces.o cf_forcep.o\
             cf_hopping.o\
             cf_induce.o\
             cf_init.o\
             cf_list.o\
             cf_main.o\
             cf_pairs.o\
             cf_pbc.o\
             cf_rdgeom.o\
             cf_shake.o\
             cf_grid.o\
             cf_rdtop.o\
             cf_uhop.o\
             pme_init.o\
             pme_fft.o\
             pme_grid.o

LIBRARY=libcafe.a

HEADER=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h
