#
# $Id: GNUmakefile,v 1.9 2004-09-30 06:46:30 marat Exp $
#

 OBJ_OPTIMIZE =  qmmm_main.o qmmm_data.o qmmm_bq_data.o mm_interface.o qmmm_geom.o  qmmm.o qmmm_input.o mm_utils.o

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
