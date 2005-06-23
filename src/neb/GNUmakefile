#
# $Id: GNUmakefile,v 1.1 2005-06-23 00:17:14 marat Exp $
#

OBJ_OPTIMIZE =  neb_input.o neb_drv.o bead_list.o

LIBRARY=libneb.a

HEADERS=

USES_BLAS= neb_drv.F

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
