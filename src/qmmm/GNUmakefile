#
# $Id: GNUmakefile,v 1.1 2000-10-23 20:10:47 d3j191 Exp $
#

 OBJ_OPTIMIZE =  qmmm_main.o

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
