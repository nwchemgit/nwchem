#
# $Id: GNUmakefile,v 1.6 2000-10-24 02:38:23 d3j191 Exp $
#

 OBJ_OPTIMIZE =  nwmd.o md_start.o md_input.o md_main.o \
                 md_finish.o \
                 prp_init.o prp_main.o 

LIBRARY=libnwmd.a

HEADERS=constants.fh md_common.fh

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
