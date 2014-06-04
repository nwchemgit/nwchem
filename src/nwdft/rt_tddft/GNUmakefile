OBJ = \
rt_tddft.o


OBJ_OPTIMIZE = 

LIBRARY = libnwdft.a

USES_BLAS =

LIB_DEFINES = -DDEBUG_PRINT

SUBDIRS = matutils rtutils input init canorg propagators closedshell openshell

##
## KAL: hack to allow symlinked directory
##
LIB_INCLUDES = -I../include -I../../ddscf
#LIB_INCLUDES = -I./headers -I$(NWCHEM_TOP)/src/nwdft/include -I$(NWCHEM_TOP)/src/ddscf

HEADERS = ./headers/rt_tddft.fh ./headers/matutils.fh

##
## KAL: hack to allow symlinked directory
##
include ../../config/makefile.h
include ../../config/makelib.h
#include $(NWCHEM_TOP)/src/config/makefile.h
#include $(NWCHEM_TOP)/src/config/makelib.h
