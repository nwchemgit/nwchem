#
# $Id: GNUmakefile,v 1.4 2004-02-12 23:14:33 marat Exp $
#

 OBJ_OPTIMIZE =  mm_interface.o qmmm.o qmmm_input.o qmmm_main.o mm_utils.o

LIBRARY=libqmmm.a

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
