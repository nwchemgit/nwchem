#
# $Id$
#

LIBRARY = libnwpython.a

OBJ = python_input.o task_python.o nwchem_wrap.o nw_inp_from_string.o 

ifndef PYTHONHOME
error:
	@echo Must define PYTHONHOME and PYTHONVERSION to build with Python
	@exit 1
endif

# Under WIN32 only need to define PYTHONHOME
ifndef PYTHONVERSION
      PYTHONVERSION=$(shell $(PYTHONHOME)/bin/python -V |cut -c 8-10)
#error1:
#	@echo Must define PYTHONHOME and PYTHONVERSION to build with Python
#	@exit 1
endif

PYMAJOR:=$(word 1, $(subst ., ,$(PYTHONVERSION)))
     GOTMINGW32= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep MINGW32|cut -c21)
ifeq ($(GOTMINGW32),1)
LIB_INCLUDES =  -I$(PYTHONHOME)/include -I$(PYTHONHOME)/Include -I$(PYTHONHOME)
else
ifeq (${PYMAJOR},3)
LIB_INCLUDES = $(shell $(PYTHONHOME)/bin/python3-config --includes)
else
LIB_INCLUDES = $(shell $(PYTHONHOME)/bin/python-config --includes)
endif
endif


include ../config/makefile.h
include ../config/makelib.h


