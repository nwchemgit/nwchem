#
# $Id: GNUmakefile,v 1.7 2005-10-03 21:18:07 edo Exp $
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
ifeq ($(TARGET),LINUX64)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
