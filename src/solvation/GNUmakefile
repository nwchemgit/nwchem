#
# $Id: GNUmakefile 25404 2014-03-30 15:43:31Z niri $
#
   OBJ_OPTIMIZE = cosmo.o \
         mnsol_mem.o\
         mnsol_interface.o\
         mnsol.o \
         solv_data.o

   LIBRARY = libsolvation.a

# If you add -DJEFF_OUTPUT to LIB_DEFINES, the dipole polarizability will be printed in the same format as the TCE.
# This is very useful (1) if you wrote scripts to parse that particular formatting and (2) if you want the anisotrophy.
#   THE LAPACK CALL ENABLED BY THIS IS BROKEN ON SOME ARCHITECTURES
#   LIB_DEFINES = -DJEFF_OUTPUT

   LIB_INCLUDES = -I../nwdft/include -I../ddscf

USES_BLAS = 

HEADERS = ../property/prop.fh ../property/atomdata.fh

include ../config/makefile.h
include ../config/makelib.h

#FDEBUG += -O0
