#
# $Id: GNUmakefile,v 1.1 2000-10-24 02:37:35 d3j191 Exp $
#

 OBJ_OPTIMIZE =  qmd_start.o qmd_main.o 

LIBRARY=libqmd.a

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
