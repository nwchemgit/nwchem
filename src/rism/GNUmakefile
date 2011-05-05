#
# $Id: GNUmakefile 19708 2010-10-29 18:04:21Z d3y133 $
#


OBJ =  task_rism.o rism_input.o

LIBRARY=librism.a

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
