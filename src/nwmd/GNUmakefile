#
# $Id: GNUmakefile,v 1.4 2000-01-06 21:20:13 d3j191 Exp $
#

 OBJ_OPTIMIZE =  nwmd.o md_start.o md_input.o md_main.o \
                 md_finish.o \
                 prp_init.o prp_main.o \
                 qmd_start.o qmd_main.o qmmm_main.o

LIBRARY=libnwmd.a

HEADERS=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
