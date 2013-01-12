#
# $Id$
#

OBJ_OPTIMIZE =  neb_input.o neb_drv.o neb_utils.o bead_list.o 

LIBRARY=liboptim.a

HEADERS=

USES_BLAS= neb_drv.F neb_utils.F bead_list.F neb_input.F

include ../../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
