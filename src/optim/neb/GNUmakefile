#
# $Id: GNUmakefile,v 1.2 2005-08-08 23:50:37 marat Exp $
#

OBJ_OPTIMIZE =  neb_input.o neb_drv.o neb_utils.o bead_list.o 

LIBRARY=liboptim.a

HEADERS=

USES_BLAS= neb_drv.F

include ../../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
