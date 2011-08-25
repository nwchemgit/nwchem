#
# $Id: GNUmakefile 19708 2010-10-29 18:04:21Z d3y133 $
#


OBJ =  task_rism.o rism_input.o rism_prepare.o db_data.o rism_thermo.o \
       rism_util.o rism_output.o

LIBRARY=librism.a

HEADERS=

USES_BLAS=

include ../config/makefile.h

#DEFINES += -DRISM_DEBUG
ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
