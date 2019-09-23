#
# $Id$
#

LIBRARY = libnwpython.a

OBJ = python_input.o task_python.o nwchem_wrap.o nw_inp_from_string.o 


# Under WIN32 only need to define PYTHONHOME
ifndef PYTHONVERSION
      GOTPYTHON3 := $(shell command -v python3 2> /dev/null)
      ifdef GOTPYTHON3
      PYTHONVERSION=$(shell python3 -V 2>&1 |cut -c 8-10)
      else
      PYTHONVERSION=$(shell python -V 2>&1 |cut -c 8-10)
      endif
endif

PYMAJOR:=$(word 1, $(subst ., ,$(PYTHONVERSION)))
#     GOTMINGW32= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep MINGW32|cut -c21)
#ifeq ($(GOTMINGW32),1)
#LIB_INCLUDES =  -I$(PYTHONHOME)/include -I$(PYTHONHOME)/Include -I$(PYTHONHOME)
#else
ifeq (${PYMAJOR},3)
LIB_INCLUDES = $(shell python3-config --includes)
else
LIB_INCLUDES = $(shell python-config --includes)
endif
#endif


include ../config/makefile.h
include ../config/makelib.h


