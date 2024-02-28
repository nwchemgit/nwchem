#
# $Id$
#

LIBRARY = libnwpython.a

OBJ = python_input.o task_python.o nwchem_wrap.o nw_inp_from_string.o 

ifndef PYTHONVERSION
   PYTHONVERSION=$(shell python -V| cut -d ' ' -f 2 |cut -d . -f 1-2)
endif

LIB_INCLUDES = $(shell python$(PYTHONVERSION)-config --includes| sed -e "s/\-iwithsysroot/\-I/g")



include ../config/makefile.h
include ../config/makelib.h


