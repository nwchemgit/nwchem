#
# $Id: GNUmakefile,v 1.3 2004-02-10 01:07:32 marat Exp $
#

 OBJ_OPTIMIZE =  mm_interface.o qmmm.o qmmm_input.o qmmm_main.o

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
