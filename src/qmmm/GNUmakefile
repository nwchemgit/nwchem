#
# $Id: GNUmakefile,v 1.2 2000-10-23 21:13:28 d3j191 Exp $
#

 OBJ_OPTIMIZE =  qmmm_input.o qmmm_main.o

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
