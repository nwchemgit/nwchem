#
# $Id: GNUmakefile,v 1.4 2000-08-26 19:57:15 d3g681 Exp $
#

LIBRARY = libpython.a

OBJ = python_input.o task_python.o nwchem_wrap.o nw_inp_from_string.o 

ifndef PYTHONHOME
error:
	@echo Must define PYTHONHOME and PYTHONVERSION to build with Python
	@exit 1
endif

# Under WIN32 only need to define PYTHONHOME
ifndef PYTHONVERSION
error1:
	@echo Must define PYTHONHOME and PYTHONVERSION to build with Python
	@exit 1
endif


LIB_INCLUDES = -I$(PYTHONHOME)/include/python$(PYTHONVERSION)


include ../config/makefile.h
include ../config/makelib.h


