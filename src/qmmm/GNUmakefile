#
# $Id: GNUmakefile,v 1.7 2004-07-21 01:29:32 marat Exp $
#

 OBJ_OPTIMIZE =  qmmm_data.o mm_interface.o qmmm_geom.o  qmmm_bq.o qmmm.o qmmm_input.o mm_utils.o

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
