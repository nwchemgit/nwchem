#
# $Id: GNUmakefile,v 1.3 1999-07-27 21:21:20 d3e129 Exp $
#

LIBRARY = libpython.a

OBJ = python_input.o task_python.o nwchem_wrap.o nw_inp_from_string.o 

ifndef PYTHONHOME
error:
	@echo Must define PYTHONHOME and PYTHONVERSION to build with Python
	@exit 1
endif
ifndef PYTHONVERSION
error1:
	@echo Must define PYTHONHOME and PYTHONVERSION to build with Python
	@exit 1
endif

LIB_INCLUDES += -I$(PYTHONHOME)/include/python$(PYTHONVERSION)


include ../config/makefile.h
include ../config/makelib.h


