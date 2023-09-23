#$Id$

OBJ_OPTIMIZE =  map_hcore.o rot_vir.o ducc.o hh.o hp.o pp.o hhhh.o pppp.o hhhp.o ppph.o hhpp.o hphp.o
LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = rot_vir.F

LIB_DEFINES = -DDEBUG_PRINT

ifdef OPTIMIZE_DUCC
LIB_DEFINES += -DOPTIMIZE_DUCC
endif

include ../../config/makefile.h
include ../../config/makelib.h

V=1

ifeq ($(FC),ifx)
#fix for hhhh.F segv crash with ifx 2023.2.0
FOPTIMIZE = -O1
endif


