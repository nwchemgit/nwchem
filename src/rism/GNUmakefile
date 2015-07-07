#
# $Id$
#


OBJ =  task_rism.o rism_input.o rism_prepare.o db_data.o \
       rism_util.o rism_output.o
# task_rism.o rism_input.o rism_prepare.o db_data.o rism_thermo.o \

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
