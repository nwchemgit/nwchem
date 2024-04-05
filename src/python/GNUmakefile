#
# $Id$
#

LIBRARY = libnwpython.a

OBJ = python_input.o task_python.o nwchem_wrap.o nw_inp_from_string.o 

ifndef PYTHONVERSION
   GOTPYTHON := $(shell command -v python 2> /dev/null)
   GOTPYTHON2 := $(shell command -v python2 2> /dev/null)
   GOTPYTHON3 := $(shell command -v python3 2> /dev/null)
   ifdef GOTPYTHON3
      PYTHONVERSION=$(shell python3 -V| cut -d ' ' -f 2 |cut -d . -f 1-2)
   else ifdef GOTPYTHON2
      PYTHONVERSION=$(shell python2 -V| cut -d ' ' -f 2 |cut -d . -f 1-2)
   else ifdef GOTPYTHON
      PYTHONVERSION=$(shell python -V| cut -d ' ' -f 2 |cut -d . -f 1-2)
   endif
endif

LIB_INCLUDES = $(shell python$(PYTHONVERSION)-config --includes| sed -e "s/\-iwithsysroot/\-I/g")



include ../config/makefile.h
include ../config/makelib.h


